/*
 * By Ray Murdorf
 * Licensed under CC-BY-NC-SA 4.0
 */

public class FreeSpline {
    static int n; // total number of points given
    static double[] supe; // size n - 3
    static double[] main; // size n - 2
    static double[] sub; // size n - 3
    static double[] constants; // size n - 2

    public static void gaussianElimination(/* double[] supe, double[] main, double[] sub, double[] constants, int n */) throws Exception {
        double mod = 0;
        // this for loop eliminates the entire subdiagonal
        // the super diagonal is not modified since every time we are subtracting 0 from it
        for(int i = 0; i < sub.length; i++) {
            mod = sub[i] / main[i];
            // the operation we are actually doing, sub[i] -= main[i] * sub[i] / main[i] always simplifies to 0
            sub[i] = 0;
            main[i + 1] -= supe[i] * mod;
            // catches if any element of the "main" diagonal becomes 0
            // if(Math.abs(main[i]) < 0.000000001) throw new Exception("Singular system");
            constants[i + 1] -= constants[i] * mod;
        }
        // checks last "main" element, in this case a_n,n
        if(Math.abs(main[main.length - 1]) < 0.000000001) throw new Exception("Singular system");
        // after this step we are left with a triangular matrix, which can be evaluated with back substitution
    }
    public static double[] backSubstitution(/* double[] supe, double[] main, double[] constants, int n */) {
        for(int i = n - 4; i >= 0; i--) {
            constants[i + 1] /= main[i + 1];
            // like above, this operation will always simplify to the same answer, in this case 1
            main[i + 1] = 1;
            constants[i] -= constants[i + 1] * supe[i];
            // again, this always simplifies to zero
            supe[i] = 0;
        }
        constants[0] /= main[0];
        main[0] = 1;
        return constants;
    }

    public static void main(String[] args) {
        double[] x = {0.0, 0.4, 1.0, 1.4, 2.0, 2.4, 3.0, 3.4, 4.0};
        double[] y = {-2.4, 0.2, 2.0, -0.5, -2.6, 0.3, 3.5, -1.2, 2.4};
        n = x.length;
        double[] a = new double[n];
        double[] b = new double[n - 1];
        double[] c = new double[n];
        double[] d = new double[n - 1];
        double[] h = new double[n - 1];
        constants = new double[n - 2];
        supe = new double[n - 3];
        sub = new double[n - 3];
        main = new double[n - 2];

        // fills the a matrix
        for(int i = 0; i < a.length; i++) {
            a[i] = y[i];
        }

        // fills h matrix
        for(int i = 0; i < h.length; i++) {
            h[i] = x[i + 1] - x[i];
        }

        // computes constants (aka k matrix)
        for(int i = 1; i < constants.length + 1; i++) {
            constants[i - 1] = ((3 * (a[i + 1] - a[i])) / h[i]) - ((3 * (a[i] - a[i - 1])) / h[i - 1]);
        }

        // computes the main, sub-, and super-diagonals for the matrix
        for(int i = 0; i < supe.length; i++) {
            supe[i] = h[i + 1];
            sub[i] = h[i + 1];
        }
        for(int i = 0; i < main.length; i++) {
            main[i] = 2 * (h[i] + h[i + 1]);
        }

        try {
            // computes c matrix
            gaussianElimination();
            double[] temp = backSubstitution();
            c[0] = 0;
            c[n - 1] = 0;
            for(int i = 1; i < c.length - 1; i++) {
                c[i] = temp[i - 1];
            }

            // computes b matrix
            for(int i = 0; i < b.length; i++) {
                b[i] = ((a[i + 1] - a[i]) / h[i]) - ((h[i] * ((2 * c[i]) + c[i + 1])) / 3);
            }

            // computes the d matrix
            for(int i = 0; i < d.length; i++) {
                d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
            }

            // prints splines
            for(int i = 0; i < x.length - 1; i++) {
                System.out.printf("S%d(x) = %.3f + %.3f(x - %.2f) + %.3f(x - %.2f)^2 + %.3f(x - %.2f)^3 %n", i, a[i], b[i], x[i], c[i], x[i], d[i], x[i]);
            }
        } catch(Exception e) {
            System.out.println(e);
        }
    }
}
