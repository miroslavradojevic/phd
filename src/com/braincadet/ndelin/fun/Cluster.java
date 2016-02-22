package com.braincadet.ndelin.fun;

import com.braincadet.ndelin.multi.X;

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

//    private static int[] clustering(ArrayList<X> x) {
//
//        int[] labels = new int[x.size()];
//        for (int i = 0; i < labels.length; i++) labels[i] = i;
//
//        for (int i = 0; i < x.size(); i++) {
//            for (int j = 0; j < x.size(); j++) {
//
//                if (i!=j) {
//
//                    double dst2 	= Math.pow(x.get(i).x-x.get(j).x, 2) + Math.pow(x.get(i).y-x.get(j).y, 2) + Math.pow(x.get(i).z-x.get(j).z, 2);
//                    double rd2 		= Math.pow(x.get(i).sig + x.get(j).sig, 2);
//
//                    if (dst2<=rd2) {  // they are neighbours
//
//                        if (labels[j]!=labels[i]) {
//
//                            int currLabel = labels[j];
//                            int newLabel  = labels[i];
//
//                            labels[j] = newLabel;
//
//                            //set all that also were currLabel to newLabel
//                            for (int k = 0; k < labels.length; k++)
//                                if (labels[k]==currLabel)
//                                    labels[k] = newLabel;
//
//                        }
//
//                    }
//
//                }
//
//            }
//
//        }
//
//        return labels; // cluster labels for each disc
//
//    }

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