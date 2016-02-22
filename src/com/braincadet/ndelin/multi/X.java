package com.braincadet.ndelin.multi;

import ij.IJ;

import java.util.ArrayList;

/**
 * Created by miroslav on 10/21/15.
 */
public class X {

    public float x;
    public float y;
    public float z;

    public float vx;
    public float vy;
    public float vz;

//    public float sig;   // gaussian cross-section standard deviation
//    public float corr;  // correlation with the gaussian cross-section template that corresponds to oriented particle (x,y,z,vx,vy,vz)
    public float tness; // tubeness at (x,y,z)
    public float w;     // density weight (phd particle weight)

    public int tag; // each particle associated with a supression tag

    public X(){}

//    public void set(float[] loc, float[] dir) {
//
//        this.x = loc[0];
//        this.y = loc[1];
//        this.z = loc[2];
//
//        this.vx = dir[0];
//        this.vy = dir[1];
//        this.vz = dir[2];
//
////        this.sig = Float.NaN;
////        this.corr = Float.NaN;
////        this.w = Float.NaN;
////        this.tness = Float.NaN;
////        this.tag = -1;
//
//    }

    public X(float[] loc, float[] dir, float w, float tness) { // float sig, float corr,

        this.x = loc[0];
        this.y = loc[1];
        this.z = loc[2];

        this.vx = dir[0];
        this.vy = dir[1];
        this.vz = dir[2];

//        this.sig = sig;
//        this.corr = corr;
        this.w = w;
        this.tness = tness;

    }

    public X(double x, double y, double z, double vx, double vy, double vz, float w, float tness) { // , float sig, float corr

        this.x = (float) x;
        this.y = (float) y;
        this.z = (float) z;

        this.vx = (float) vx;
        this.vy = (float) vy;
        this.vz = (float) vz;

//        this.sig = sig;
//        this.corr = corr;
        this.w = w;
        this.tness = tness;

    }

    public X(float x, float y, float z, float vx, float vy, float vz, float w, float tness) { // , float sig, float corr

        this.x = x;
        this.y = y;
        this.z = z;

        this.vx = vx;
        this.vy = vy;
        this.vz = vz;

//        this.sig = sig;
//        this.corr = corr;
        this.w = w;
        this.tness = tness;

    }

    public X(X Xinit) {

        this.x = Xinit.x;
        this.y = Xinit.y;
        this.z = Xinit.z;

        this.vx = Xinit.vx;
        this.vy = Xinit.vy;
        this.vz = Xinit.vz;

//        this.sig = Xinit.sig;
//        this.corr = Xinit.corr;
        this.w = Xinit.w;
        this.tness = Xinit.tness;

    }

    public void set(X Xset) {

        this.x = Xset.x;
        this.y = Xset.y;
        this.z = Xset.z;

        this.vx = Xset.vx;
        this.vy = Xset.vy;
        this.vz = Xset.vz;

//        this.sig = Xset.sig;
//        this.corr = Xset.corr;
        this.w = Xset.w;
        this.tness = Xset.tness;

    }

    public static void print(X x) {
        String s = "";
        s += IJ.d2s(x.x, 2) + ", " + IJ.d2s(x.y, 2) + ", " + IJ.d2s(x.z, 2) + " | " +
                IJ.d2s(x.vx, 2) + ", " + IJ.d2s(x.vy, 2) + ", " + IJ.d2s(x.vz,2) + " | " +
//                IJ.d2s(x.corr,2) + ", " + IJ.d2s(x.sig,2) +
                IJ.d2s(x.tness,2) + " | " +
                IJ.d2s(x.w,2);
        IJ.log(s);
    }

    public static void sortByTness(ArrayList<X> x) {

        for (int i = 0; i < x.size()-1; i++) {
            for (int j = i+1; j < x.size(); j++) {
                if (x.get(j).tness>x.get(i).tness) { // descend.

                    // swap the values
                    X temp 	= new X(x.get(i));
                    x.get(i).set(x.get(j));
                    x.get(j).set(temp);

                }
            }
        }

    }

}
