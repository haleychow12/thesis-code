#include <SoftwareSerial.h>

#define MESSAGE_HEADER 'A'
#define MESSAGE_FOOTER '\n'

#define AT_DIRECTION_0 -3000
#define AT_DIRECTION_1 -1500
#define AT_DIRECTION_2 0
#define AT_DIRECTION_3 1500
#define AT_DIRECTION_4 3000

#define AT_DIRECTION_0_PIN 9
#define AT_DIRECTION_1_PIN 8
#define AT_DIRECTION_2_PIN 7
#define AT_DIRECTION_3_PIN 6
#define AT_DIRECTION_4_PIN 5

#define AT_DISTANCE_NUM_SHIFT_CHIPS 2
#define AT_DISTANCE_DATA_WIDTH (AT_DISTANCE_NUM_SHIFT_CHIPS * 8)
#define AT_DISTANCE_PULSE_WIDTH_USEC 5
#define AT_DISTANCE_POLL_DELAY_MSEC 1

#define AT_DISTANCE_PLOAD_PIN 10
#define AT_DISTANCE_CLOCKENABLE_PIN 14
#define AT_DISTANCE_DATA_PIN 15
#define AT_DISTANCE_CLOCK_PIN 16

#define AT_DISTANCE_DIGIT_MASK 0xFE
#define AT_DISTANCE_DECIMAL_MASK 0x01

#define AT_DISTANCE_DIGIT_0 0xFC // 11111100
#define AT_DISTANCE_DIGIT_1 0x60 // 01100000
#define AT_DISTANCE_DIGIT_2 0xDA // 11011010
#define AT_DISTANCE_DIGIT_3 0xF2 // 11110010
#define AT_DISTANCE_DIGIT_4 0x66 // 01100110
#define AT_DISTANCE_DIGIT_5 0xB6 // 10110110
#define AT_DISTANCE_DIGIT_6 0xBE // 10111110
#define AT_DISTANCE_DIGIT_7 0xE0 // 11100000
#define AT_DISTANCE_DIGIT_8 0xFE // 11111110
#define AT_DISTANCE_DIGIT_9 0xF6 // 11110110
#define AT_DISTANCE_DIGIT_S 0xB6 // 10110110
#define AT_DISTANCE_DIGIT_E 0x9E // 10011110

int fieldDirection;
unsigned int consecutiveDirection;

int distanceToVictim;
unsigned int consecutiveDistance;

/* setup: Configure 12C interface. */

void setup() {
  /*EDIT: Serial 1 -> Serial*/
  Serial.begin(115200);
  while (!Serial);

  //setup AT direction pins
  pinMode(AT_DIRECTION_0, INPUT);
  pinMode(AT_DIRECTION_1, INPUT);
  pinMode(AT_DIRECTION_2, INPUT);
  pinMode(AT_DIRECTION_3, INPUT);
  pinMode(AT_DIRECTION_4, INPUT);

  //setup AT distance pins
  pinMode(AT_DISTANCE_PLOAD_PIN, OUTPUT);
  pinMode(AT_DISTANCE_CLOCKENABLE_PIN, OUTPUT);
  pinMode(AT_DISTANCE_DATA_PIN, INPUT);
  pinMode(AT_DISTANCE_CLOCK_PIN, OUTPUT);

  //setup shift registers
  digitalWrite(AT_DISTANCE_CLOCK_PIN, LOW);
  digitalWrite(AT_DISTANCE_PLOAD_PIN, HIGH);

  //initialzie magnetic field variables
  fieldDirection = 0;
  distanceToVictim = 32767;

}

void loop() {
  // update magnetic field measurements
  measureDirection();
  measureDistance();

  //send data over serial
  /*EDIT: Serial1 -> Serial*/
  Serial.print(MESSAGE_HEADER);
  Serial.print(",");
  Serial.print(fieldDirection, DEC);
  Serial.print(",");
  Serial.print(distanceToVictim, DEC);
  Serial.print(",");
  Serial.print(MESSAGE_FOOTER);

}

/*measure Direction: measures the direction of the field line
 in centidegress and updates fieldDriection accordingly*/
void measureDirection() {
  //TRUE AND HIGH are defined to be 1. FALSE and LOW are defined to be 0

  //read AT's direction LEDs
  boolean directionLED[5];
  directionLED[0] = digitalRead(AT_DIRECTION_0_PIN);
  directionLED[1] = digitalRead(AT_DIRECTION_1_PIN);
  directionLED[2] = digitalRead(AT_DIRECTION_2_PIN);
  directionLED[3] = digitalRead(AT_DIRECTION_3_PIN);
  directionLED[4] = digitalRead(AT_DIRECTION_4_PIN);

  //does this reading match the previous reading?
  static boolean directionLEDPrev[5];
  boolean consecutive = directionLEDPrev[0] == directionLED[0] && 
  directionLEDPrev[1] == directionLED[1] && 
  directionLEDPrev[2] == directionLED[2] && 
  directionLEDPrev[3] == directionLED[3] && 
  directionLEDPrev[4] == directionLED[4];

  //count consecutive readings
  if (consecutive) {
    consecutiveDirection++;
  }
  else {
    consecutiveDirection = 0;
    directionLEDPrev[0] = directionLED[0];
    directionLEDPrev[1] = directionLED[1];
    directionLEDPrev[2] = directionLED[2];
    directionLEDPrev[3] = directionLED[3];
    directionLEDPrev[4] = directionLED[4];
  }

  //only interpret direction when it has stabilized
  if (consecutiveDirection < 100) {
    return;
  }

  //calculate number of LEDs turned on 
  int numOn = 0;
  for (int i = 0; i < 5; i++) {
    numOn += directionLED[i];
  }

  //return immediately if the direction cannot be determined
  if (numOn <= 0) {
    return;
  }

  //calculate sum of angles indicated by direction LEDs
  int sumAngles = directionLED[0] * AT_DIRECTION_0 +
  directionLED[1] * AT_DIRECTION_1 +
  directionLED[2] * AT_DIRECTION_2 +
  directionLED[3] * AT_DIRECTION_3 +
  directionLED[4] * AT_DIRECTION_4;

  //calculate average angle indicated by direction LEDs
  int avgAngle = sumAngles / numOn;

  //update fieldDirection 
  fieldDirection = avgAngle;
}

/*measure Distance: measures the distance along the field line to
the victim (in cm) and updates fieldDistance accordingly*/

void measureDistance() {
  
  //read shift register to capture display segment values
  unsigned int bytesVal = readShiftRegisters();

  if (bytesVal == 0) {
    consecutiveDistance = 0;
    return;
  }

  //count consecutive readings
  static unsigned int bytesValPrev;
  if (bytesValPrev == bytesVal) {
    consecutiveDistance++;
  }
  else {
    consecutiveDistance = 0;
    bytesValPrev = bytesVal;
  }

  //only interpret distance display when it has stabilized
  if (consecutiveDistance < 150) {
    return;
  }

  //read left digit 
  byte leftDigitRaw = (bytesVal >> 8) & AT_DISTANCE_DIGIT_MASK;
  boolean leftDecimal = (bytesVal >> 8) & AT_DISTANCE_DECIMAL_MASK;

  //read right digit
  byte rightDigitRaw = bytesVal & AT_DISTANCE_DIGIT_MASK;
  boolean rightDecimal = bytesVal & AT_DISTANCE_DECIMAL_MASK;

  //maintain previous distance if display shows "SE"
  int leftS = (leftDigitRaw == AT_DISTANCE_DIGIT_S);
  int rightE = (rightDigitRaw == AT_DISTANCE_DIGIT_E);

  if (leftS && rightE) {
    return;
  }

  //interpret left digit value
  byte leftDigit;
  switch (leftDigitRaw) {
    case AT_DISTANCE_DIGIT_0:
      leftDigit = 0;
      break;
    case AT_DISTANCE_DIGIT_1:
      leftDigit = 1;
      break;
    case AT_DISTANCE_DIGIT_2:
      leftDigit = 2;
      break;
    case AT_DISTANCE_DIGIT_3:
      leftDigit = 3;
      break;
    case AT_DISTANCE_DIGIT_4:
      leftDigit = 4;
      break;
    case AT_DISTANCE_DIGIT_5:
      leftDigit = 5;
      break;
    case AT_DISTANCE_DIGIT_6:
      leftDigit = 6;
      break;
    case AT_DISTANCE_DIGIT_7:
      leftDigit = 7;
      break;
    case AT_DISTANCE_DIGIT_8:
      leftDigit = 8;
      break;
    case AT_DISTANCE_DIGIT_9:
      leftDigit = 9;
      break;
    default:
      leftDigit = 0;
      break;
  }

  byte rightDigit;
  switch (rightDigitRaw) {
    case AT_DISTANCE_DIGIT_0:
      rightDigit = 0;
      break;

    case AT_DISTANCE_DIGIT_1:
      rightDigit = 1;
      break;

    case AT_DISTANCE_DIGIT_2:
      rightDigit = 2;
      break;

    case AT_DISTANCE_DIGIT_3:
      rightDigit = 3;
      break;

    case AT_DISTANCE_DIGIT_4:
      rightDigit = 4;
      break;

    case AT_DISTANCE_DIGIT_5:
      rightDigit = 5;
      break;

    case AT_DISTANCE_DIGIT_6:
      rightDigit = 6;
      break;

    case AT_DISTANCE_DIGIT_7:
      rightDigit = 7;
      break;

    case AT_DISTANCE_DIGIT_8:
      rightDigit = 8;
      break;

    case AT_DISTANCE_DIGIT_9:
      rightDigit = 9;
      break;

    default:
      rightDigit = 0;
      break;
  }

  //interpret decimal placement
  if (leftDecimal) {
    distanceToVictim = (int) (100* (leftDigit + rightDigit/10.0));
  }
  else {
    distanceToVictim = 100 * (10 * leftDigit + rightDigit);
  }

}

/*read shift registers: triggers parallel load on at distance num shift chips
shift register. then shifts in each loaded value and returns the state of the 
shift registers as an unsigned int. */

unsigned int readShiftRegisters(){
  int bitVal;
  unsigned int bytesVal = 0;

  //trigger a parallel load to latch the state of the data lines
  digitalWrite(AT_DISTANCE_CLOCKENABLE_PIN, HIGH);
  digitalWrite(AT_DISTANCE_PLOAD_PIN, LOW);
  delayMicroseconds(AT_DISTANCE_PULSE_WIDTH_USEC);
  digitalWrite(AT_DISTANCE_PLOAD_PIN, HIGH);
  digitalWrite(AT_DISTANCE_CLOCKENABLE_PIN, LOW);

  //read each latched value on the serial output line
  for (int i = 0; i < AT_DISTANCE_DATA_WIDTH; i++){
    bitVal = digitalRead(AT_DISTANCE_DATA_PIN);
    bytesVal |= (bitVal << i);

    //pulse clock to shift-in next bit
    digitalWrite(AT_DISTANCE_CLOCKENABLE_PIN, HIGH);
    delayMicroseconds(AT_DISTANCE_PULSE_WIDTH_USEC);
    digitalWrite(AT_DISTANCE_CLOCKENABLE_PIN, LOW);
  }

  return bytesVal;
}

