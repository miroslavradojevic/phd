package com.braincadet.phd.fun;

import com.braincadet.phd.multi.X;
import ij.IJ;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

/**
 * Created by miroslav on 10/13/15.
 */
public class Tools {

    public static void createAndCleanDir(String dirpath) {
        File f1 = new File(dirpath);
        if (!f1.exists()) {f1.mkdirs();}
        else { // clean it, reset with each exec
            File[] files = f1.listFiles();
            for (File file : files)
                if (!file.delete()) System.out.println("Failed to delete " + file);
        }
//        IJ.log(dirpath);
    }

    public static void createDir(String dirpath) {
        // create directory without cleaning it up
        File f1 = new File(dirpath);
        if (!f1.exists()) {f1.mkdirs();}
//        IJ.log("createDir " + dirpath);
    }

    public static void cleanfile(String filepath) {

        try {
            PrintWriter logWriter = new PrintWriter(filepath);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

    }

    public static float[][] calcdirs180(boolean is2D, int nrdirs) { // sample uniformly on semisphere

        float[][] vxyz = new float[nrdirs][3]; // set of unit directions expressed with 3d vectors

        double h_k, theta_k, phi_k, phi_k_1 = 0;

        for (int k = 0; k < nrdirs; k++) { // generate Ndir directions

            if (is2D) {
                float ang1 = 0f + k * (3.14f/(float)nrdirs);
                vxyz[k][0] = (float) Math.cos(ang1);
                vxyz[k][1] = (float) Math.sin(ang1);
                vxyz[k][2] =  0;
            }
            else {

                h_k = 1 - 1 * ((double)k/(nrdirs-1)); // 1 : 0 defines angular range
                theta_k = Math.acos(h_k);

                if(k==0 || k==(nrdirs-1)) {
                    phi_k   = 0;
                    phi_k_1 = 0;
                }
                else {
                    phi_k = phi_k_1 + 3.6 / ( Math.sqrt(nrdirs) * Math.sqrt(1 - h_k * h_k));
                    phi_k_1 = phi_k;
                }

                vxyz[k][0] = (float) (Math.sin(theta_k) * Math.cos(phi_k));
                vxyz[k][1] = (float) (Math.sin(theta_k) * Math.sin(phi_k));
                vxyz[k][2] = (float)  Math.cos(theta_k);

            }

        }

        return vxyz;

    }

    public static double[][] calcdirs360(boolean is2D, int nrdirs) { // sample uniformly on sphere

        double[][] vxyz = new double[nrdirs][3]; // set of unit directions expressed with 3d vectors

        double h_k, theta_k, phi_k, phi_k_1 = 0;

        for (int k = 0; k < nrdirs; k++) { // generate Ndir directions

            if (is2D) {
                float ang1 = 0f + k * ((2*3.14f)/(float)nrdirs);
                vxyz[k][0] = (float) Math.cos(ang1);
                vxyz[k][1] = (float) Math.sin(ang1);
                vxyz[k][2] =  0;
            }
            else {

                h_k = 1 - 2 * ((double)k/(nrdirs-1)); // 1 : -1 defines angular range
                theta_k = Math.acos(h_k);

                if(k==0 || k==(nrdirs-1)) {
                    phi_k   = 0;
                    phi_k_1 = 0;
                }
                else {
                    phi_k = phi_k_1 + 3.6 / ( Math.sqrt(nrdirs) * Math.sqrt(1 - h_k * h_k));
                    phi_k_1 = phi_k;
                }

                vxyz[k][0] = (float) (Math.sin(theta_k) * Math.cos(phi_k));
                vxyz[k][1] = (float) (Math.sin(theta_k) * Math.sin(phi_k));
                vxyz[k][2] = (float)  Math.cos(theta_k);

            }

        }

        return vxyz;

    }

    public static void rotation_matrix(float a1, float a2, float b1, float b2, float[][] R) {
        // from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        // assume a(a1,a2) and b(b1,b2) are unit vectors

        // v is the norm of the cross product of (a1,a2) and (b1,b2)
        float v = a1*b2-a2*b1; // direction is z but it does not exist here
        float dot_prod = a1*b1+a2*b2;
        float cross_prod_norm_2 = v*v;

        if (cross_prod_norm_2<=0.00001) {

            if (Math.abs(dot_prod - 1)<Math.abs(dot_prod + 1)) {
                // identity mapping (a and b are aligned)
                R[0][0] = 1;
                R[0][1] = 0;
                R[1][0] = 0;
                R[1][1] = 1;
            }
            else {
                // inversion (a and b are opposite)
                R[0][0] = -1;
                R[0][1] = 0;
                R[1][0] = 0;
                R[1][1] = -1;
            }

        }
        else {
            // cross product is not small - safe to calculate
            R[0][0] = dot_prod;
            R[0][1] = -v;
            R[1][0] = v;
            R[1][1] = dot_prod;
        }

    }

    public static void rotation_apply(float[][] R, float v1, float v2, float[] out12) {
        out12[0] = R[0][0]*v1 + R[0][1]*v2;
        out12[1] = R[1][0]*v1 + R[1][1]*v2;
    }

    public static void rotation_matrix(float a1, float a2, float a3, float b1, float b2, float b3, float[][] R) {
        // from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        // assume a(a1,a2,a3) and b(b1,b2,b3) are unit vectors

        // v is cross product of (a1, a2, a3) and (b1, b2, b3)
        float v1 = a2*b3 - b2*a3;
        float v2 = -(a1*b3-b1*a3);
        float v3 = a1*b2-b1*a2;

        // cross product is zero - symmetric operations
        float cross_prod_norm_2     = v1*v1+v2*v2+v3*v3;
        float dot_prod              = a1*b1+a2*b2+a3*b3;

        if (cross_prod_norm_2<=0.00001) { // cross product is negligable

            if (Math.abs(dot_prod - 1)<Math.abs(dot_prod + 1)) {
                // identity mapping (a and b are aligned)
                R[0][0] = 1;
                R[0][1] = 0;
                R[0][2] = 0;

                R[1][0] = 0;
                R[1][1] = 1;
                R[1][2] = 0;

                R[2][0] = 0;
                R[2][1] = 0;
                R[2][2] = 1;
            }
            else {
                // inversion (a and b are opposite)
                R[0][0] = -1;
                R[0][1] = 0;
                R[0][2] = 0;

                R[1][0] = 0;
                R[1][1] = -1;
                R[1][2] = 0;

                R[2][0] = 0;
                R[2][1] = 0;
                R[2][2] = -1;
            }
        }
        else {
            // cross product is not negligable - safe to calculate
            float tt = (1-(a1*b1+a2*b2+a3*b3))/(v1*v1+v2*v2+v3*v3);
            R[0][0] = 1 + 0     + tt * (-v3*v3-v2*v2);
            R[0][1] = 0 + (-v3) + tt * (v1*v2);
            R[0][2] = 0 + (v2)  + tt * (-v1*v3);

            R[1][0] = 0 + (v3)  + tt * (v1*v2);
            R[1][1] = 1 + 0     + tt * (-v3*v3-v1*v1);
            R[1][2] = 0 + (-v1) + tt * (v2*v3);

            R[2][0] = 0 + (-v2) + tt * (v1*v3);
            R[2][1] = 0 + (v1)  + tt * (v2*v3);
            R[2][2] = 1 + 0     + tt * (-v2*v2-v1*v1);
        }

    }

    public static void rotation_apply(float[][] R, float v1, float v2, float v3, float[] out123) {

        out123[0] = R[0][0]*v1 + R[0][1]*v2 + R[0][2]*v3;
        out123[1] = R[1][0]*v1 + R[1][1]*v2 + R[1][2]*v3;
        out123[2] = R[2][0]*v1 + R[2][1]*v2 + R[2][2]*v3;

    }

    public static float interp(float atX, float atY, float atZ, float[] img, int width, int height, int length) {

        int x1 = (int) atX;
        int x2 = x1 + 1;
        float x_frac = atX - x1;

        int y1 = (int) atY;
        int y2 = y1 + 1;
        float y_frac = atY - y1;

        if (length==1) { // atZ is not necessary

            boolean isIn2D = x1>=0 && x2<width && y1>=0 && y2<height;
            if(!isIn2D) return 0;//Float.NaN;

            // take neighbourhood 2d
            float I11_1 = img[y1*width+x1];
            float I12_1 = img[y1*width+x2];
            float I21_1 = img[y2*width+x1];
            float I22_1 = img[y2*width+x2];

            return (1-y_frac) * ((1-x_frac)*I11_1 + x_frac*I12_1) + (y_frac) * ((1-x_frac)*I21_1 + x_frac*I22_1);

        }
        else {

            int z1 = (int) atZ;
            int z2 = z1 + 1;
            float z_frac = atZ - z1;

            boolean isIn3D = y1>=0 && y2<height && x1>=0 && x2<width && z1>=0 && z2<length;
            if (!isIn3D) return 0;//Float.NaN;

            // take neighbourhood 3d
            float I11_1 = img[z1*width*height+y1*width+x1];
            float I12_1 = img[z1*width*height+y1*width+x2];
            float I21_1 = img[z1*width*height+y2*width+x1];
            float I22_1 = img[z1*width*height+y2*width+x2];

            float I11_2 = img[z2*width*height+y1*width+x1];
            float I12_2 = img[z2*width*height+y1*width+x2];
            float I21_2 = img[z2*width*height+y2*width+x1];
            float I22_2 = img[z2*width*height+y2*width+x2];

            return (1-z_frac)  *
                    (  (1-y_frac) * ((1-x_frac)*I11_1 + x_frac*I12_1) + (y_frac) * ((1-x_frac)*I21_1 + x_frac*I22_1) )   +
                    z_frac      *
                            (  (1-y_frac) * ((1-x_frac)*I11_2 + x_frac*I12_2) + (y_frac) * ((1-x_frac)*I21_2 + x_frac*I22_2) );

        }


    }

    public static float interp(X x, float[] img, int imgw, int imgh, int imgl) {
        return interp(x.x, x.y, x.z, img, imgw, imgh, imgl);
    }

    public static int[] descending(ArrayList<Integer> a) {

        // prepare array with indexes first
        int[] idx = new int[a.size()];
        for (int i=0; i<idx.length; i++) idx[i] = i;

        for (int i = 0; i < a.size()-1; i++) {
            for (int j = i+1; j < a.size(); j++) {
                if (a.get(j)>a.get(i)) { // desc.
                    int temp 	= a.get(i);
                    a.set(i, a.get(j));
                    a.set(j, temp);

                    int temp_idx 	= idx[i];
                    idx[i] 			= idx[j];
                    idx[j]			= temp_idx;
                }
            }
        }

        return idx;

    }

    public static int[] asc(ArrayList<Float> a) {

        int[] idx = new int[a.size()];
        for (int i=0; i<idx.length; i++) idx[i] = i;

        for (int i = 0; i < a.size()-1; i++) {
            for (int j = i+1; j < a.size(); j++) {
                if (a.get(j)<a.get(i)) { // asc.

                    float temp 	= a.get(i);
                    a.set(i, a.get(j));
                    a.set(j, temp);

                    int temp_idx 	= idx[i];
                    idx[i] 			= idx[j];
                    idx[j]			= temp_idx;
                }
            }
        }

        return idx;

    }

    public static int[] desc(ArrayList<Float> a) {

        int[] idx = new int[a.size()];
        for (int i=0; i<idx.length; i++) idx[i] = i;

        for (int i = 0; i < a.size()-1; i++) {
            for (int j = i+1; j < a.size(); j++) {
                if (a.get(j)>a.get(i)) { // desc.

                    float temp 	= a.get(i);
                    a.set(i, a.get(j));
                    a.set(j, temp);

                    int temp_idx 	= idx[i];
                    idx[i] 			= idx[j];
                    idx[j]			= temp_idx;
                }
            }
        }

        return idx;

    }

    static void pick(float[] cdf, float[][] Pxyz, float[] Qxyz) {

        // Pxyz.length particles will be used to pick 1 using cdf pty distribution
        // systematic resampling algorithm

        int i = 0;

        float uj = new Random().nextFloat(); // uniform [0,1)

        while (uj > cdf[i]) {
            i++;
        }

        Qxyz[0] = Pxyz[i][0];
        Qxyz[1] = Pxyz[i][1];
        Qxyz[2] = Pxyz[i][2];

    }

    static void pick(float[] csw, float[][] Pxyz, float[][] Qxyz) {

        // Pxyz.length particles will be used to resample Qxyz.length ones from csw (pdf)
        if (Qxyz.length>Pxyz.length) {
            Qxyz = null;
            return;
        }

        if (Qxyz.length==Pxyz.length) {
            // take all particles

            return;
        }

        // select Qxyz.length samples
        int np = Qxyz.length;
        int i = 0;

        float u1 = (1f/(float)np) * new Random().nextFloat(); // uniform [0,1/np)

        for (int j = 0; j < np; j++) {

            float uj = u1 + j/(float)np;

            while (uj > csw[i]) {
                i++;
            }

            Qxyz[j][0] = Pxyz[i][0];
            Qxyz[j][1] = Pxyz[i][1];
            Qxyz[j][2] = Pxyz[i][2];

        }

    }

    static float quantile(float[] a, int ratioNum, int ratioDen) {

        int n = a.length;
        int i, j, l, m, k;
        double x;

        if (ratioNum>=ratioDen) k = n-1;
        else k = (int) Math.floor(n * ((float) ratioNum / (float) ratioDen));
//    else if ((ratioNum*n) % ratioDen == 0) k = ((ratioNum*n)/ratioDen)-1;
//    else k = (ratioNum*n)/ratioDen;

        l=0 ; m=n-1;
        while (l < m) {
            x=a[k];
            i = l;
            j = m;
            do {
                while (a[i] < x) i++ ;
                while (x < a[j]) j-- ;
                if (i <= j) {
                    float temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                    i++ ; j-- ;
                }
            } while (i <= j);
            if (j < k) l = i;
            if (k < i) m = j;
        }

        return a[k];

    }

    public static void quicksort(float[] main, int[] index) {
        quicksort(main, index, 0, index.length - 1);
    }

    // quicksort a[left] to a[right]
    public static void quicksort(float[] a, int[] index, int left, int right) {
        if (right <= left) return;
        int i = partition(a, index, left, right);
        quicksort(a, index, left, i-1);
        quicksort(a, index, i+1, right);
    }

    // partition a[left] to a[right], assumes left < right
    private static int partition(float[] a, int[] index, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (less(a[++i], a[right]))      // find item on left to swap
                ;                               // a[right] acts as sentinel
            while (less(a[right], a[--j]))      // find item on right to swap
                if (j == left) break;           // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(a, index, i, j);               // swap two elements into place
        }
        exch(a, index, i, right);               // swap with partition element
        return i;
    }

    // is x < y ?
    private static boolean less(float x, float y) {
        return (x < y);
    }

    // exchange a[i] and a[j]
    private static void exch(float[] a, int[] index, int i, int j) {
        float swap = a[i];
        a[i] = a[j];
        a[j] = swap;
        int b = index[i];
        index[i] = index[j];
        index[j] = b;
    }

    public static final void normalize(ArrayList<Float> in) { // will normalize with respect to the max.

        float curr_max = in.get(0);

        for (int i = 1; i < in.size(); i++) {
            if(in.get(i)>curr_max) curr_max = in.get(i);
        }

        if (curr_max>2*Float.MIN_VALUE) {
            for (int i=0; i<in.size(); i++) {
                in.set(i, in.get(i)/curr_max);
            }
        }


    }

    /**
     * background median estimation ( Wirth's algorithm )
     * Title: Algorithms + data structures = programs
     * Publisher: Englewood Cliffs: Prentice-Hall, 1976
     * Physical description: 366 p.
     * Series: Prentice-Hall Series in Automatic Computation
     */
    public static double median_Wirth(float[] a) {
        int n = a.length;
        int i, j, l, m, k;
        double x;
        if (n % 2 == 0) k = (n/2)-1;
        else k = (n/2);
        l=0 ; m=n-1 ;
        while (l < m)
        {
            x=a[k] ;
            i = l ;
            j = m ;
            do
            {
                while (a[i] < x) i++ ;
                while (x < a[j]) j-- ;
                if (i <= j) {
                    float temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                    i++ ; j-- ;
                }
            } while (i <= j) ;
            if (j < k) l = i ;
            if (k < i) m = j ;
        }
        return a[k] ;
    }

}
