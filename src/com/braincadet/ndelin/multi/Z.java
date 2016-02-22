package com.braincadet.ndelin.multi;

/**
 * Created by miroslav on 10/25/15.
 */
public class Z {

    public float x;
    public float y;
    public float z;
//    public float    sig; // cross-section gaussian standard deviation will correspond to the size of the observation
    public int      count; // how many out of ro*np are >correlation threshold, will be used to determine if it is clutter

    public Z(){}

    public Z(float x, float y, float z, int count) { // float sig,

        this.x = x;
        this.y = y;
        this.z = z;
//        this.sig = sig;
        this.count = count;

    }

    public Z(Z zin) {

        this.x = zin.x;
        this.y = zin.y;
        this.z = zin.z;
//        this.sig = zin.sig;
        this.count = zin.count;

    }

//    public Z(Xk Xin) {
//
//        this.x      = Xin.x;
//        this.y      = Xin.y;
//        this.z      = Xin.z;
//        this.sig    = Xin.sig;
//        this.count = 0;
//
//    }

//    public Z(Xk xin, int count) {
//
//        this.x = xin.x;
//        this.y = xin.y;
//        this.z = xin.z;
//        this.sig = xin.sig;
//        this.count = count;
//
//    }

}
