import java.util.*;
public class Point{
    //Cartesian stuff will be done in this class
    //magnetic moment 
    private static final double MOMENT = 0.02714336053;
    private static double incX, incY;

    private double x, y, e;
    private int d;
    public Point(double xPoint, double yPoint){
        this.x = xPoint;
        this.y = yPoint;
        this.d = -1;
    }

    private void setError(double e, int degree){
        this.e = e;
        this.d = degree;
    }

    private static void setIncrementers(double xincr, double yincr){
        incX = xincr;
        incY = yincr;
    }


    //rotates x and y around 0,0 by theta radians, returns point with new x and y
    private static Point rotate(double x, double y, double T){
        return new Point(x*Math.cos(T) - y*Math.sin(T), x*Math.sin(T) + y*Math.cos(T));

    }

    //rotates x,y about rotate_point_x, rotate_point_y by theta radians
    //returns point with new x and y
    private static Point rotate_point(double x, double y, double rotate_point_x,
                             double rotate_point_y, double T){
        Point prime = rotate(x-rotate_point_x, y-rotate_point_y, T);

        return new Point(prime.x + rotate_point_x, prime.y + rotate_point_y);
    }

    public static Point step(Point tloc, Point myloc, double theta, double stepsize){
        //to do.
        double lx = myloc.x;
        double ly = myloc.y;
        double[] results = new double[2];
        field(tloc, myloc, theta, results);
        double fx = results[0];
        double fy = results[1];
        double c = stepsize/(Math.sqrt(fx*fx + fy*fy));

        return new Point(lx+c*fx, ly+c*fy);
    }

    //calculates and returns the distance between two points
    private static double distance(Point a, Point b){
        return (Math.sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)));
    }

    //returns an array with the value of the field from a transmitter at the
    //location tloc with orientation theta degrees, at point myloc. x value 
    //of the field located in results[0], y value in results[1]
    private static void field(Point tloc, Point myloc, double theta, double[] results){
        theta = Math.toRadians(theta);
        Point prime = rotate_point(myloc.x, myloc.y, tloc.x, tloc.y, -theta);
        
        double[] r = {prime.x - tloc.x, prime.y - tloc.y};
        double rnorm = Math.sqrt(r[0]*r[0] + r[1]*r[1]);

        double f1 = 1 / (4 * Math.PI);

        double[] f2 = {0.0, 0.0};
        double[] f3 = {0.0, 0.0};

        f2[0] = 3 * r[0] * (MOMENT*r[1]) / (Math.pow(rnorm, 5));
        f2[1] = 3 * r[1] * (MOMENT*r[1]) / (Math.pow(rnorm, 5));
        f3[1] = MOMENT / (rnorm*rnorm*rnorm); 

        results[0] = f1 * (f2[0] - f3[0]) * 1E6;
        results[1] = f1 * (f2[1] - f3[1]) * 1E6;

    }

    //creates a Points 2-D array that contains the locations of each of the points in the grid
    private static Point[][] fillPointsList(double minx, double miny, int xPoints, int yPoints){
        Point[][] pointList = new Point[xPoints][yPoints];

        double xincr = -1*(minx/(xPoints/2 + 1));
        double yincr = -1*(miny/(yPoints/2 + 1));

        setIncrementers(xincr, yincr);

        double x = minx;
        double y = miny;
        
        for (int i = 0; i < xPoints; i++){
            for (int j = 0; j < yPoints; j++){
                pointList[i][j] = new Point(x, y);
                //System.out.println(String.format("Point (%.4f, %.4f)", x, y));
                y += yincr;
            }
            x += xincr;
            y = miny;
        }

        return pointList;
    }

    private static double snorm(double error){
        double NORM = 20.0;
        if (error > NORM)
            return 1;
        return error/NORM;
    }

    private static double dnorm(double error){
        double NORM = 10.0;
        if (error > NORM)
            return 1;
        return error/NORM;
    }

    private static Point findRandomNeighbor(int x, int y, Point[][] pointsList, int[] nextIndex){
        int counter = 0;
        Point[] neighbors = new Point[8];
        int[] xVals = new int[8];
        int[] yVals = new int[8];

        int lenx = pointsList.length;
        int leny = pointsList[0].length;

        if (x > 0){
            neighbors[counter] = pointsList[x-1][y];
            xVals[counter] = x-1;
            yVals[counter] = y;
            counter++;
            if (y > 0){
                neighbors[counter] = pointsList[x-1][y-1];
                xVals[counter] = x-1;
                yVals[counter] = y-1;
                counter++;
            }
            if (y+1 < leny){
                neighbors[counter] = pointsList[x-1][y+1];
                xVals[counter] = x-1;
                yVals[counter] = y+1;
                counter++;
            }
        }
        if (x+1 < lenx){
            neighbors[counter] = pointsList[x+1][y];
            xVals[counter] = x+1;
            yVals[counter] = y;
            counter++;
            if (y > 0){
                neighbors[counter] = pointsList[x+1][y-1];
                xVals[counter] = x+1;
                yVals[counter] = y-1;
                counter++;
            }
            if (y+1 < leny){
                neighbors[counter] = pointsList[x+1][y+1];
                xVals[counter] = x+1;
                yVals[counter] = y+1;
                counter++;
            }
        }
        if (y > 0){
            neighbors[counter] = pointsList[x][y-1];
            xVals[counter] = x;
            yVals[counter] = y-1;
            counter++;
        }
        if (y+1 < leny){
            neighbors[counter] = pointsList[x][y+1];
            xVals[counter] = x;
            yVals[counter] = y+1;
            counter++;
        }

        int random = (int) (Math.random()*counter);
        nextIndex[0] = xVals[random];
        nextIndex[1] = yVals[random];
        return neighbors[random];
    }

    private static void storeGuess(ArrayList<Point> bestGuess, double error, Point p, int degree){
        int bestGuessLength = 4;
        int count = bestGuess.size();
        Point point;

        if (count < bestGuessLength || error < bestGuess.get(count-1).e){
            point = new Point(p.x, p.y);
            point.setError(error, degree);
            bestGuess.add(point);

            // Sorting
            Collections.sort(bestGuess, new Comparator<Point>() {
                @Override
                public int compare(Point a, Point b)
                {
                    if (a.e > b.e)
                        return 1;
                    if (a.e == b.e)
                        return 0;
                    else 
                        return -1;
                }
            });

            if (count+1 > bestGuessLength)
                bestGuess.remove(count);
        }
        //System.out.println(String.format("Source: %.4f, %.4f, Error: %.4f, Degree: %d", p.x,p.y,error, degree));

    }

    private static Point findAvg(ArrayList<Point> bestGuess){
        double avgx = 0;
        double avgy = 0;

        double errorThreshold = .75; //really need to change this, like 3?
        Point ret;

        for (Point p: bestGuess){
            avgx += p.x;
            avgy += p.y;

            if (p.e > errorThreshold){
                System.out.println("Error is too large");
                return null;
            }
        }

        double xguess = avgx/bestGuess.size();
        double yguess = avgy/bestGuess.size();

        return new Point(xguess, yguess);
    }

    //calculate the error with a transmitter at testPoint oriented in the direction theta
    //degrees
    private static double calcError(double theta, Point testPoint, Point[] searchList, 
        double[] dirList, double[] rList){
        int NUMVALS = 7; //only test the last 7 values
        int SLOPE_MULTIPLIER = 16; //multiplier for the direciton difference

        double slopeError = 0;
        double distError = 0;

        double[] testResults = new double[2];

        int length = searchList.length;
        if (length < NUMVALS)
            NUMVALS = length;


        //calculating error on only the last NUMVALS measurements
        for (int i = 0; i < NUMVALS; i++){
            int k = length - (i+1);
            field(testPoint, searchList[k], theta, testResults);

            double testr = distance(testPoint, searchList[k]);
            double testslope = testResults[1]/testResults[0];
            
            slopeError += (dirList[k] - testslope)*(dirList[k] - testslope);
            distError += (rList[k] - testr)*(rList[k] - testr);
        }
        return Math.sqrt(SLOPE_MULTIPLIER*snorm(slopeError) + dnorm(distError));
    }

    public static Point annealingAlgorithm(Point[] searchList, double[] dirList, double[] rList){
        ArrayList<Point> bestGuess = new ArrayList<Point>();
        double minx = -40;
        double miny = -20;
        int xPoints = 1000;
        int yPoints = 1000;
        Point[][] pointsList = fillPointsList(minx, miny, xPoints, yPoints);

        int interval = 2;
        int samples = 500;
        Point startPoint = null;

        //set annealing constants
        double alpha = .9;
        int jmax = 5000;
        double errormax = .01;


        //run this for some discretized amt of rotations
        for (int degree = 0; degree < 180; degree += interval){
            double temp = 9000;
            double errormin = Integer.MAX_VALUE;
            double olderror = 0;
            int[] nextIndex = new int[2];
            int xIndex = 0; int yIndex = 0;
            //find lowest value in (sample#) of samples
            for (int count = 0; count < samples; count++){
                int testxIndex = (int) (Math.random()*xPoints);
                int testyIndex = (int) (Math.random()*yPoints);
                Point testPoint = pointsList[testxIndex][testyIndex];

                olderror = calcError(degree, testPoint, searchList, dirList, rList);
                pointsList[testxIndex][testyIndex].setError(olderror, degree);

                //switch out values if olderror is less than current min
                if (olderror < errormin){
                    errormin = olderror;
                    startPoint = testPoint;
                    xIndex = testxIndex;
                    yIndex = testyIndex;
                }
            }

            //System.out.println(String.format("Source: %.4f, %.4f, Error: %.4f, Degree: %d", 
            //startPoint.x, startPoint.y, startPoint.e, startPoint.d));

            int j = 1;
            while (j <= jmax && olderror > errormax){
                //System.out.println("point:" + Double.toString(startPoint.x) + "," + Double.toString(startPoint.y));
                //System.out.println("still point:" + Double.toString(pointsList[xIndex][yIndex].x) + ", " +Double.toString(pointsList[xIndex][yIndex].y));
                double nexterror = 0;

                //get a random neighbor
                Point nextPoint = findRandomNeighbor(xIndex, yIndex, pointsList, nextIndex);           

                if (nextPoint.d == degree)
                    nexterror = nextPoint.e;
                else {
                    //calculate error
                    nexterror = calcError(degree, nextPoint, searchList, dirList, rList);
                    nextPoint.setError(nexterror, degree);
                    storeGuess(bestGuess, nexterror, nextPoint, degree);
                }

                //System.out.println(String.format("Source: %.4f, %.4f, Error: %.4f, Degree: %d", 
                //nextPoint.x, nextPoint.y, nextPoint.e, nextPoint.d));

                //simulated annealing
                double delta = nexterror - olderror;
                if (delta < 0){
                    startPoint = nextPoint;
                    olderror = nexterror;
                    xIndex = nextIndex[0];
                    yIndex = nextIndex[1];
                }
                else{

                    double p = Math.exp(-delta/temp);
                    //System.out.println(Double.toString(p));
                    if (Math.random() < p){
                        startPoint = nextPoint;
                        olderror = nexterror;
                        xIndex = nextIndex[0];
                        yIndex = nextIndex[1];
                    }
                }
                temp = temp*alpha;
                //System.out.println("temp: " + Double.toString(temp));
                j+=1;
            }
        }
        for (int i = 0; i < bestGuess.size(); i++){
            System.out.println(String.format("Guess: %.4f, %.4f, Error: %.4f, Degree: %d", 
            bestGuess.get(i).x, bestGuess.get(i).y, bestGuess.get(i).e, bestGuess.get(i).d));

        }
        return findAvg(bestGuess);

    }

    public static void main(String args[]){
        //System.out.println(Double.toString(snorm(10)));
        //System.out.println(Double.toString(dnorm(7)));

        Point[] searchList = {new Point(15.0, -7.5),
            new Point(14.9666518861, -7.3027998395),
            new Point(14.9316034054, -7.10589477608),
            new Point(14.8948381785, -6.90930287792),
            new Point(14.8563825119, -6.71303479161),
            new Point(14.8162254863, -6.51710773044),
            new Point(14.7743985392, -6.32153107823)};

        double[] dirList = {-5.91338272621, -5.61807699865,-5.34962183837,
            -5.1044730899,-4.87968382989, -4.67278364413, -4.4816860498};

        double[] rList = {31.146239984, 31.055684673, 30.9646141949, 30.8730280491,
            30.7809257092, 30.688306624, 30.5951702165};

        Point source = annealingAlgorithm(searchList, dirList, rList);
        System.out.println(String.format("Source: %.4f, %.4f", 
            source.x, source.y, source.e));


        // double e = calcError(80, new Point(5, 2), searchList, dirList, rList);
        // System.out.println(Double.toString(e));
        // Point [][] pointsList = fillPointsList(-20, -10, 500, 500);
        // Point startPoint = new Point(19.92, 9.56);
        // int x = (int) ((startPoint.x+20)/incX);
        // int y = (int) ((startPoint.y+10)/incY);
        // System.out.println(Integer.toString(x));
        // System.out.println(Integer.toString(y));

        // System.out.println(Double.toString(pointsList[x][y].x) + ","
        //           + Double.toString(pointsList[x][y].y));

        // findRandomNeighbor(x,y, pointsList);

        // ArrayList <Point> bestGuess = new ArrayList<Point>();
        // storeGuess(bestGuess, .1, new Point(15, 12), 30);
        // storeGuess(bestGuess, .5, new Point(10, 12), 30);
        // storeGuess(bestGuess, .7, new Point(13, 12), 30);
        // storeGuess(bestGuess, .2, new Point(11, 12), 30);
        // storeGuess(bestGuess, .3, new Point(14, 12), 30);
        // storeGuess(bestGuess, .03, new Point(17, 12), 30);

        // for (int i = 0; i < bestGuess.size(); i++){
        //     System.out.println(Double.toString(bestGuess.get(i).x) + "," + Double.toString(bestGuess.get(i).y));
        // }

        // for (int i = 0; i < 500; i++){
        //     for (int j = 0; j < 500; j++){
        //         System.out.println(Double.toString(pointsList[i][j].x) + ","
        //          + Double.toString(pointsList[i][j].y));
        //     }

        // }

        //double[] results = new double[2];

        //field(new Point(0,1), new Point(5,2), 45, results);
        //double[] results = new double[2];
        // Point one = rotate(3.0,6.0, (Math.PI/4));
        // System.out.println(Double.toString(one.x) + "," + Double.toString(one.y));

        // Point two = rotate_point(4.0,6.0,1.0,1.0,(Math.PI/3));
        // System.out.println(Double.toString(two.x) + "," + Double.toString(two.y));
    }
}