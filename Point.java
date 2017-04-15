import java.util.*;
public class Point {
    //Cartesian stuff will be done in this class
    //magnetic moment 
    private static final double MOMENT = 0.02714336053;

    private double x, y;
    private double e;
    int d;
    public Point(double xPoint, double yPoint){
        this.x = xPoint;
        this.y = yPoint;
        this.d = -1;
    }

    private void setError(double e, int degree){
        this.e = e;
        this.d = degree;
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

        double xincr = -1*(minx/(xPoints/2));
        double yincr = -1*(miny/(yPoints/2));
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
            int k = length - NUMVALS;
            field(testPoint, searchList[k], theta, testResults);

            double testr = distance(testPoint, searchList[k]);
            double testslope = testResults[1]/testResults[0];

            slopeError += (dirList[k] - testslope)*(dirList[k] - testslope);
            distError += (rList[k] - testr)*(rList[k] - testr);
        }
        return Math.sqrt(SLOPE_MULTIPLIER*snorm(slopeError) + dnorm(distError));
    }

    public static void annealingAlgorithm(Point[] searchList, double[] dirList, double[] rList){
        ArrayList<Point> bestGuess = new ArrayList<Point>();
        double minx = -40;
        double miny = -20;
        int xPoints = 500;
        int yPoints = 500;
        Point[][] pointsList = fillPointsList(minx, miny, xPoints, yPoints);

        int interval = 2;
        int samples = 50;

        //set annealing constants
        double alpha = .9;
        int jmax = 5000;
        double errormax = .0001;

        for (int degree = 0; degree < 180; degree += interval){
            double T = 9000;
            double errormin = Integer.MAX_VALUE;
            for (int count = 0; count < samples; count++){
                int testxIndex = (int) (Math.random()*xPoints);
                int testyIndex = (int) (Math.random()*yPoints);
                Point testPoint = pointsList[testxIndex][testyIndex];

                double error = calcError(degree, testPoint, searchList, dirList, rList);
            }
        }


    }

    public static void main(String args[]){
        //System.out.println(Double.toString(snorm(10)));
        //System.out.println(Double.toString(dnorm(7)));

        Point[] searchList = {new Point(14.9666518861, -7.3027998395),
            new Point(14.9316034054, -7.10589477608),
            new Point(14.8948381785, -6.90930287792),
            new Point(14.8563825119, -6.71303479161),
            new Point(14.8162254863, -6.51710773044),
            new Point(14.7743985392, -6.32153107823),
            new Point(14.7308435468, -6.12633127639)};

        double[] dirList = {-5.91338272621, -5.61807699865,-5.34962183837,
            -5.1044730899,-4.87968382989, -4.67278364413, -4.4816860498};

        double[] rList = {31.055684673, 30.9646141949, 30.8730280491,
            30.7809257092, 30.688306624, 30.5951702165, 30.5015158843};


        double e = calcError(66.12, new Point(-14.7067, 1.8595), searchList, dirList, rList);
        System.out.println(Double.toString(e));
        // Point [][] pointsList = fillPointsList(-20, -10, 500, 500);

        // for (int i = 0; i < 500; i++){
        //     System.out.println(Double.toString(pointsList[i][0].x));
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