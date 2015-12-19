package com.braincadet.ndelin.fun;

import com.braincadet.ndelin.multi.X;
import ij.IJ;

import java.util.ArrayList;

/**
 * Created by miroslav on 11/8/15.
 */
public class Cluster {

    public static int   MAXITER     = Integer.MAX_VALUE;
    public static float EPSILON2    = 0.000001f;

//    public static void meanShiftX2Z(ArrayList<X> xin, float kdist, int minCount, int maxNobj, ArrayList<Z> zout) {
//        float[][] convxyz = meanShift(xin, MAXITER, EPSILON2, kdist);
//        int[] lab = clustering(convxyz, 1f);
//        extractZ(lab, xin, minCount, maxNobj, zout);
//    }

    public static float[][] meanShift(ArrayList<X> xin, float dist) {

        float[][] conv = new float[xin.size()][3];
        for (int i = 0; i < xin.size(); i++) {
            conv[i][0] = xin.get(i).x;
            conv[i][1] = xin.get(i).y;
            conv[i][2] = xin.get(i).z;
        }

        float[] new_v = new float[3]; // auxiliary variable for iteration (to avoid allocation inside the loop)

        for (int i = 0; i < conv.length; i++) {

            int iter = 0;
            double d2;

            do {
                runOne(conv[i], new_v, xin, dist);

                d2 = Math.pow(new_v[0]-conv[i][0],2) + Math.pow(new_v[1]-conv[i][1],2) + Math.pow(new_v[2]-conv[i][2],2);

                conv[i][0] = new_v[0];
                conv[i][1] = new_v[1];
                conv[i][2] = new_v[2];

                iter++;
            }
            while (iter < MAXITER && d2 > EPSILON2);

        }

        return conv;

    }

    public static float[][] meanShift(ArrayList<X> xin, float[] xw, float dist) {

        float[][] conv = new float[xin.size()][3];

        for (int i = 0; i < xin.size(); i++) {
            conv[i][0] = xin.get(i).x;
            conv[i][1] = xin.get(i).y;
            conv[i][2] = xin.get(i).z;
        }

        float[] new_v = new float[3];

        for (int i = 0; i < conv.length; i++) {

            int iter = 0;
            double d2;

            do {
                runOne(conv[i], new_v, xin, xw, dist);

                d2 = Math.pow(new_v[0]-conv[i][0],2) + Math.pow(new_v[1]-conv[i][1],2) + Math.pow(new_v[2]-conv[i][2],2);

                conv[i][0] = new_v[0];
                conv[i][1] = new_v[1];
                conv[i][2] = new_v[2];

                iter++;
            }
            while (iter < MAXITER && d2 > EPSILON2);

        }

        return conv;

    }

    private static void runOne(float[] curr_v, float[] new_v, ArrayList<X> xin, float[] xw, float dist){

        float sum = 0;
        new_v[0] = 0;
        new_v[1] = 0;
        new_v[2] = 0;

        for (int l = 0; l < xin.size(); l++) {

            float x2 = curr_v[0]-xin.get(l).x; x2 = x2*x2;
            float y2 = curr_v[1]-xin.get(l).y; y2 = y2*y2;
            float z2 = curr_v[2]-xin.get(l).z; z2 = z2*z2;

            if (x2+y2+z2 <= dist*dist) {

                sum += xw[l];

                new_v[0] += xin.get(l).x * xw[l];
                new_v[1] += xin.get(l).y * xw[l];
                new_v[2] += xin.get(l).z * xw[l];

            }

        }

        if (sum>0) {

            new_v[0] /= sum;
            new_v[1] /= sum;
            new_v[2] /= sum;

        }

    }

    private static void runOne(float[] curr_v, float[] new_v, ArrayList<X> xin, float dist) {

        float sum 	= 0;
        new_v[0] 	= 0;
        new_v[1] 	= 0;
        new_v[2]    = 0;

        for (int l = 0; l < xin.size(); l++) {

            float x2 = curr_v[0]-xin.get(l).x; x2 = x2*x2;
            float y2 = curr_v[1]-xin.get(l).y; y2 = y2*y2;
            float z2 = curr_v[2]-xin.get(l).z; z2 = z2*z2;

            if (x2+y2+z2 <= dist*dist) {

                sum += 1;

                new_v[0] += xin.get(l).x;
                new_v[1] += xin.get(l).y;
                new_v[2] += xin.get(l).z;

            }

        }

        if (sum>0) {

            new_v[0] /= sum;
            new_v[1] /= sum;
            new_v[2] /= sum;

        }
//        else {
//            IJ.log("happened, " + dist);
//            new_v[0] = curr_v[0];
//            new_v[1] = curr_v[1];
//            new_v[2] = curr_v[2];
//        }

    }

    public static int[] clustering(float[][] values, float dist) {
        // indxs represent indexes of values that need to be clustered according to their values read before
        // dists are the distances (idxs.length * idxs.length),
        // threshold_dists is the distance limit
        // output is list of unique labels

        int[] labels = new int[values.length];
        for (int i = 0; i < labels.length; i++) labels[i] = i; // each object gets its unique label

        float[][] dists = new float[values.length][values.length]; // not really efficient but calculate it
        for (int i = 0; i < values.length; i++) {
            for (int j = i; j < values.length; j++) {
                if(i==j) {
                    dists[i][j] = 0;
                }
                else {
                    float x2 = values[i][0]-values[j][0]; x2 = x2*x2;
                    float y2 = values[i][1]-values[j][1]; y2 = y2*y2;
                    float z2 = values[i][2]-values[j][2]; z2 = z2*z2;

                    dists[i][j] = x2+y2+z2;
                    dists[j][i] = dists[i][j];
                }
            }
        }

        for (int i = 0; i < values.length; i++) {

            // one versus the rest
            for (int j = 0; j < values.length; j++) {

                if (i != j) {

                    if (dists[i][j]<=dist*dist) {

                        if (labels[j] != labels[i]) {

                            int currLabel = labels[j];
                            int newLabel  = labels[i];

                            labels[j] = newLabel;

                            //set all that also were currLabel to newLabel
                            for (int k = 0; k < labels.length; k++)
                                if (labels[k]==currLabel)
                                    labels[k] = newLabel;

                        }

                    }

                }

            }

        }

        return labels;

    }

    private static int[] clustering(ArrayList<X> x) {

        int[] labels = new int[x.size()];
        for (int i = 0; i < labels.length; i++) labels[i] = i;

        for (int i = 0; i < x.size(); i++) {
            for (int j = 0; j < x.size(); j++) {

                if (i!=j) {

                    double dst2 	= Math.pow(x.get(i).x-x.get(j).x, 2) + Math.pow(x.get(i).y-x.get(j).y, 2) + Math.pow(x.get(i).z-x.get(j).z, 2);
                    double rd2 		= Math.pow(x.get(i).sig + x.get(j).sig, 2);

                    if (dst2<=rd2) {  // they are neighbours

                        if (labels[j]!=labels[i]) {

                            int currLabel = labels[j];
                            int newLabel  = labels[i];

                            labels[j] = newLabel;

                            //set all that also were currLabel to newLabel
                            for (int k = 0; k < labels.length; k++)
                                if (labels[k]==currLabel)
                                    labels[k] = newLabel;

                        }

                    }

                }

            }

        }

        return labels; // cluster labels for each disc

    }

    public static void estimate(int[] labels, ArrayList<X> xin, int count_min, int nclust_max, ArrayList<X> xout) {

        boolean[] checked = new boolean[labels.length];

        ArrayList<X> xclust = new ArrayList<X>();
        ArrayList<Integer> cnt = new ArrayList<Integer>();

        for (int i = 0; i < labels.length; i++) {

            if (!checked[i]) {

                float wsum = xin.get(i).w;
                float c_x = xin.get(i).x * xin.get(i).w;
                float c_y = xin.get(i).y * xin.get(i).w;
                float c_z = xin.get(i).z * xin.get(i).w;

                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < labels.length; j++) {
                    if (!checked[j] && labels[j]==labels[i]) {

                            wsum += xin.get(j).w;
                            c_x += xin.get(j).x * xin.get(j).w;
                            c_y += xin.get(j).y * xin.get(j).w;
                            c_z += xin.get(j).z * xin.get(j).w;

                            count++;
                            checked[j] = true;

                    }
                }

                if (count >= count_min) {
                    if (wsum>0) {
                        xclust.add(new X(c_x/wsum, c_y/wsum, c_z/wsum, Float.NaN, Float.NaN, Float.NaN, 1f, Float.NaN, 1f, Float.NaN));
                        cnt.add(count);
                    }
                    else {
                        IJ.log("there was cluster with enough counts but zero wsum!");
                    }
                }

            }
        }

        int[] desc_idx = Tools.descending(cnt); // cnt list will be modified too

        xout.clear();
//        xcount.clear();

//        IJ.log("estimating "+Math.min(nclust_max,xclust.size())+" clusters:");

        for (int ii=0; ii<Math.min(nclust_max,xclust.size()); ii++) {

//            IJ.log("cluster" + ii + ": count=" + cnt.get(ii) + "   x=" + xclust.get(desc_idx[ii]).x);
            xout.add(new X(xclust.get(desc_idx[ii])));
//            xcount.add(cnt.get(ii));

        }

    }

//    public static float reweight(int[] labels, float[][] conv_xyz, ArrayList<X> xin, ArrayList<X> xclust, float K, int clutter_count) {
//
//        boolean[] checked = new boolean[labels.length];
//
//        xclust.clear();
//        float mass = 0;
//
//        for (int i = 0; i < labels.length; i++) {
//
//            if (!checked[i]) {
//
//                // add another cluster
//                ArrayList<Integer> idxs = new ArrayList<Integer>();
//
//                float wsum =                          xin.get(i).w; // used weighted mean with weights that are
//                float xe =      conv_xyz[i][0]      * xin.get(i).w; // xin.get(i).x
//                float ye =      conv_xyz[i][1]      * xin.get(i).w; // xin.get(i).y
//                float ze =      conv_xyz[i][2]      * xin.get(i).w; // xin.get(i).z
//                float vxe =     xin.get(i).vx       * xin.get(i).w;
//                float vye =     xin.get(i).vy       * xin.get(i).w;
//                float vze =     xin.get(i).vz       * xin.get(i).w;
//                float tnesse =  xin.get(i).tness    * xin.get(i).w;
////                int count = 1;
//                idxs.add(i);
//                checked[i] = true;
//
//                // check the rest
//                for (int j = i+1; j < labels.length; j++) {
//                    if (!checked[j] && labels[j]==labels[i]) {
//                        wsum   +=                      xin.get(j).w;
//                        xe     +=   conv_xyz[j][0]   * xin.get(j).w;
//                        ye     +=   conv_xyz[j][1]   * xin.get(j).w;
//                        ze     +=   conv_xyz[j][2]   * xin.get(j).w;
//                        vxe    +=   xin.get(j).vx    * xin.get(j).w;
//                        vye    +=   xin.get(j).vy    * xin.get(j).w;
//                        vze    +=   xin.get(j).vz    * xin.get(j).w;
//                        tnesse +=   xin.get(j).tness * xin.get(j).w;
////                        count++;
//                        idxs.add(j);
//                        checked[j] = true;
//                    }
//                }
//
////                IJ.log(labels[i] + " : "+wsum + " " + Float.MIN_VALUE);
//
//                // X(x, y, z, vx, vy, vz, sig, corr, w, tness)
//                if (wsum>Float.MIN_VALUE) {
//                    xclust.add(new X(xe / wsum, ye / wsum, ze / wsum, vxe / wsum, vye / wsum, vze / wsum, 1, 1, wsum / idxs.size(), tnesse / idxs.size()));
//                }
//
//                // cluster has N elements
//                int N = idxs.size();
//                float cltr = MultiTT.clutter(N, K, clutter_count);
//                float wup = 1f/(N*(1+cltr));
//
//                mass += 1f/(1f+cltr);
//
//                // assign it to xin
//                for (int j = 0; j < N; j++) {
//                    xin.get(idxs.get(j)).w = wup; // reweight PHD particles from xin
//                }
//
//            }
//        }
//
//        return mass;
//
//    }

    public static void extract(int[] labels, float[][] conv, int count_min, int nclust_max, ArrayList<X> xout, ArrayList<Integer> xcount) {

        boolean[] checked = new boolean[labels.length];

        ArrayList<X> xclust = new ArrayList<X>();
        ArrayList<Integer> cnt = new ArrayList<Integer>();

        for (int i = 0; i < labels.length; i++) {

            if (!checked[i]) {

                float cx2  = conv[i][0];
                float cy2  = conv[i][1];
                float cz2  = conv[i][2];

                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < labels.length; j++) {
                    if (!checked[j] && labels[j]==labels[i]) {

                        cx2 += conv[j][0];
                        cy2 += conv[j][1];
                        cz2 += conv[j][2];

                        count++;
                        checked[j] = true;

                    }
                }

                if (count >= count_min) {
                    xclust.add(new X(cx2/count, cy2/count, cz2/count, Float.NaN, Float.NaN, Float.NaN, 1f, Float.NaN, 1f, Float.NaN));
                    cnt.add(count);
                }
            }
        }

        int[] desc_idx = Tools.descending(cnt); // cnt will be modified too

        xout.clear();
        xcount.clear();

//        IJ.log("extracting "+Math.min(nclust_max,xclust.size())+" clusters:");

        for (int ii=0; ii<Math.min(nclust_max,xclust.size()); ii++) {

//            IJ.log("cluster" + ii + ": count=" + cnt.get(ii));
            xout.add(new X(xclust.get(desc_idx[ii])));
            xcount.add(cnt.get(ii));

        }

    }

//    private static void extractZ(int[] labels, ArrayList<X> xin, int min_count, int maxNclust, ArrayList<Z> zout) {
//
//        boolean[] checked = new boolean[labels.length];
//        ArrayList<Integer>  cnt         = new ArrayList<Integer>();
//        ArrayList<Z>        zclust      = new ArrayList<Z>();
//
//        for (int i = 0; i < labels.length; i++) {
//
//            if (!checked[i]) {
//
//                float c_x = xin.get(i).x;
//                float c_y = xin.get(i).y;
//                float c_z = xin.get(i).z;
//                float c_r = xin.get(i).sig;
//
//                int count = 1;
//                checked[i] = true;
//
//                // check the rest
//                for (int j = i+1; j < labels.length; j++) {
//                    if (!checked[j]) {
//                        if (labels[j]==labels[i]) {
//
//                            c_x += xin.get(j).x;
//                            c_y += xin.get(j).y;
//                            c_z += xin.get(j).z;
//                            c_r += xin.get(j).sig;
//
//                            count++;
//                            checked[j] = true;
//
//                        }
//                    }
//                }
//
//                if (count >= min_count) {
//                    zclust.add(new Z(c_x/count, c_y/count, c_z/count, c_r/count, count));
//                    cnt.add(count);
//                }
//
//            }
//        }
//
//        int[] desc_idx = Tools.descending(cnt); // cnt will be modified too
//
////        String s = "";
////        for (int i = 0; i < cnt.size(); i++)
////            s += " " + cnt.get(i) + " ";
////        IJ.log( ""+cnt.size() + " clusters :   " + s);
//
//        zout.clear();
//
//        for (int ii=0; ii<Math.min(maxNclust,zclust.size()); ii++)
//            zout.add(new Z(zclust.get(desc_idx[ii])));
//
//    }

}