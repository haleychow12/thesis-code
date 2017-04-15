import java.util.*;
public class Point {
    //Cartesian stuff will be done in this class
    //magnetic moment 
    private static final double MOMENT = 0.02714336053;

    private double x, y;
    private double e, d;
    public Point(double xPoint, double yPoint){
        x = xPoint;
        y = yPoint;
        d = -1;
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

    public static void annealingAlgorithm(Point[] searchList, double[] dirList, double[] rList){
        ArrayList<Point> bestGuess = new ArrayList<Point>();
        double minx = -40;
        double miny = -20;
        int xPoints = 500;
        int yPoints = 500;
        Point[][] pointsList = fillPointsList(minx, miny, xPoints, yPoints);

        int interval = 2;
        int samples = 50;

        for (int degree = 0; degree < 180; degree += interval){
            double errormin = Integer.MAX_VALUE;
            for (int count = 0; count < samples; count++){
                int testxIndex = (int) (Math.random()*xPoints);
                int testyIndex = (int) (Math.random()*yPoints);
                Point testPoint = pointsList[testxIndex][testyIndex];

                //double error = calcError(testPoint, searchList, dirList, rList);
            }
        }


    }

    public static void main(String args[]){

        // Point [][] pointsList = fillPointsList(-20, -10, 500, 500);

        // for (int i = 0; i < 500; i++){
        //     System.out.println(Double.toString(pointsList[i][0].x));
        // }

         double[] results = new double[2];

         field(new Point(0,1), new Point(5,2), 45, results);
        //double[] results = new double[2];
        // Point one = rotate(3.0,6.0, (Math.PI/4));
        // System.out.println(Double.toString(one.x) + "," + Double.toString(one.y));

        // Point two = rotate_point(4.0,6.0,1.0,1.0,(Math.PI/3));
        // System.out.println(Double.toString(two.x) + "," + Double.toString(two.y));
    }
}