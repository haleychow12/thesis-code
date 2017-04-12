public class Point {

    //Cartesian stuff will be done in this class
    double x, y;
    public Point(double xPoint, double yPoint){
        x = xPoint;
        y = yPoint;
    }

    //rotates x and y around 0,0 by theta radians
    //places x' in results[0] and y' in results[1]
    public static void rotate(double x, double y, double T, double[] results){
        results[0] = x*Math.cos(T) - y*Math.sin(T);
        results[1] = x*Math.sin(T) + y*Math.cos(T);
    }

    //rotates x,y about rotate_point_x, rotate_point_y by theta radians
    //places x' in results[0] and y' in results[1]
    public static void rotate_point(double x, double y, double rotate_point_x,
                             double rotate_point_y, double T, double[] results){
        rotate(x-rotate_point_x, y-rotate_point_y, T, results);

        results[0] = results[0] + rotate_point_x;
        results[1] = results[1] + rotate_point_y;
    }


    public static void crawlingAlgorithm(){

    }

    public static void main(String args[]){
        double[] results = new double[2];
        rotate(5.0,6.0, Math.PI, results);
        System.out.println(Double.toString(results[0]) + "," + Double.toString(results[1]));

        rotate_point(5.0,6.0,1.0,1.0,(Math.PI/2), results);
        System.out.println(Double.toString(results[0]) + "," + Double.toString(results[1]));
    }
}