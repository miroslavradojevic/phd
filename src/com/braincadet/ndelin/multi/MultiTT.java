package com.braincadet.ndelin.multi;

/*
Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

import com.braincadet.ndelin.fun.Cluster;
import com.braincadet.ndelin.fun.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class MultiTT {

    Random rndgen = new Random();

    static float gcsstd_min = 1f;           // minimal gaussian cross section std.
    static float sig2rad = 2.5f;            // how gaussian cross section relates to the tube radius
    static float gcsstd_step = .5f;         // gaussian cross section standard dev. - step
    static int   ndirs2d = 20;              // 1/2 circle, radius=step
    static int   ndirs3d = 40;              // 1/2 sphere, radius=step

    boolean     is2d;                       // image 2d or 3d
    float       gcsstd_max;                 //
    int         gcsstd_nr;                  //
    float[]     gcsstd;                     //

    int         nobjstart = -1;             // number of objects to start with
    int         ro = -1;                    // samples per particle when doing the prediction for each particle
    int         ni = -1;                    // predictions per sample

    int         diam = -1;                  // tube diameter in pixels
    int         step = -1;                  // radius for sampling new particles

    float       kappa = Float.NaN;          // von Mises distribution kappa = 0.5, 1, 2
    float       pS = Float.NaN;             // probability
    float       pD = Float.NaN;             // detection probability
    float       cluttertness = Float.NaN;   // 0-background, 1-tubularity level at which init objects are
    float       kclutt = Float.NaN;         // drop in PHD measure of the clutter
    double      tnessinit = Double.NaN;      // tubularity level at initial objects

    float[][]   vxyzUniform;                // 2d, 3d directions

    int         U, W, V, U2, W2, V2;        // ...
    float[]     img_vals;                   // ...
    float[][]   tt;                         // templates
    float[]     tta;                        // average for each template

    // for orthogonals in znccX
    float       ux, uy, uz;
    float       wx, wy, wz;

    public Stepper mm; // prediction model (offsets and weights used for importance sampling)

    // phd filtering variables
    public ArrayList<X> Xk;                     // multi object phd particles
    public ArrayList<X> XPk;                    // persistent object particle
    public ArrayList<X> Zk;                     // measurements
    public ArrayList<Integer> Zk_cnt;           // measurement counts
    public ArrayList<Float> XPk_cws;            // cummulative weight sum for the predicted particles
    public ArrayList<X> Yk;                     // estimations (from Xk) at iteration
    public ArrayList<ArrayList<X>> Y;           // accumulated iterations Y.get(i) ~ ArrayList<X>

//    public static int ITER_LIMIT = 100;         // cannot log after this limit
    public static int OBJECT_LIMIT = 500;       //

    private X[]         Xpred;   // auxilliary storage for predicted particles
    public int          npcles;  // number of particles approximating the probability density
    public float        phdmass; //

    public float[][]    g;
    public float[]      Cz;

    public MultiTT(boolean is2d, int nobjstart, int ro, int ni, int diam, int step, float kappa, float pS, float pD, float cluttertness, float kclutt) {

        this.is2d = is2d;
        this.nobjstart = nobjstart;
        this.ro = ro;
        this.ni = ni;
        this.diam = diam;
        this.step = step;
        this.kappa = kappa;
        this.pS = pS;
        this.pD = pD;
        this.cluttertness = cluttertness;
        this.kclutt = kclutt;

        mm = new Stepper(this.step, this.is2d, this.kappa, this.ro);

        Xk = new ArrayList<X>();            //
        XPk = new ArrayList<X>();           // predicted particles
        Zk = new ArrayList<X>();            // thresholded particles
        Zk_cnt = new ArrayList<Integer>();  //
        XPk_cws = new ArrayList<Float>();   //
        Yk = new ArrayList<X>();            //
        Y = new ArrayList<ArrayList<X>>();  //

        Xpred = new X[2];
        Xpred[0] = new X();
        Xpred[1] = new X();

        // g, Cz will be allocated each iteration
        vxyzUniform = Tools.calcdirs180(this.is2d, (this.is2d ? ndirs2d : ndirs3d)); // form the directions

        // templates formation (templates used to calculate the likelihood)
        // template size U*W*V, W is added for 3d (it is 1 in 2d), V is longitudinal
        U2 = this.diam/2;
        U = 2*U2 + 1; // 1st orthogonal

        if (this.is2d) {
            W2 = 0;
            W = 1;
        }
        else {
            W2 = this.diam/2;
            W = 2*W2 + 1;
        }

        // gcsstd define
        gcsstd_max = (U2/sig2rad);

        // refer V to gcsstd_max
//        V2 = this.d/2;//(this.is2d)?d/2:1;
//        V = 2*V2 + 1;

//        V2 = 1;      // + (int) (Math.ceil(gcsstd_max)/2);//(this.is2d)?d/2:1;
        V2 = this.diam/4;
        V = 2*V2 + 1;

        int cnt = 0;
        for (float sg = gcsstd_min; sg <= gcsstd_max; sg+=gcsstd_step) cnt++;
        gcsstd_nr = cnt;

        gcsstd = new float[gcsstd_nr];
        cnt = 0;
        for (float sg = gcsstd_min; sg <= gcsstd_max; sg+=gcsstd_step)
            gcsstd[cnt++] = sg;

        img_vals = new float[U*W*V];

        // single, tta templates define
        tt = new float[gcsstd_nr][]; // cross-section template
        for (int i = 0; i < gcsstd_nr; ++i) {
            tt[i] = new float[U*W*V];
        }

        tta = new float[gcsstd_nr];

        // calculate the templates
        for (int i = 0; i < gcsstd_nr; ++i) {

            float ag = 0;

            // indexing is different in 2d and 3d
            if (this.is2d) {

                for (int vv = -V2; vv <= V2; ++vv) {
                    for (int uu = -U2; uu <= U2; ++uu) {
                        //for (int ww = -W2; ww <= W2; ++ww) { // W2=0, W=1

                        float value = (float) Math.exp(-(uu * uu) / (2 * (gcsstd[i] * gcsstd[i]))); // +pow(ww,2)

                        tt[i][(vv+V2)*U+(uu+U2)] = value;
                        ag += value;

                        //}
                    }
                }

            }
            else {

                for (int vv = -V2; vv <= V2; ++vv) {
                    for (int uu = -U2; uu <= U2; ++uu) {
                        for (int ww = -W2; ww <= W2; ++ww) {

                            float value = (float) Math.exp(-(uu*uu+ww*ww)/(2*(gcsstd[i]*gcsstd[i])));

                            tt[i][(vv+V2)*U*W+(ww+W2)*U+(uu+U2)] = value;
                            ag += value;

                        }
                    }
                }

            }

            tta[i] = ag / (U * W * V);

        }

    }

    void exporttemplates(String outdir) {

        for (int i = 0; i < tt.length; i++) {

            String name = outdir+ File.separator+"template,g="+IJ.d2s(gcsstd[i],1)+".tif";

            if (this.is2d) { //
                ImageStack isout = new ImageStack(U,V); // W=1
                FloatProcessor fpout = new FloatProcessor(U, V, tt[i].clone());
                isout.addSlice(fpout);
                ImagePlus ipout = new ImagePlus("tt2d,g="+IJ.d2s(gcsstd[i],1), isout);
                IJ.saveAsTiff(ipout, name);
                IJ.log("saving:\t"+name);

            }
            else {

                float[][] ttres = new float[V][U*W];

                for (int j = 0; j < tt[i].length; j++) {
                    ttres[j/(U*W)][j%(U*W)] = tt[i][j];
                }

                ImageStack isout = new ImageStack(U,W);

                for (int j = 0; j < ttres.length; j++) {
                    FloatProcessor fpout = new FloatProcessor(U, W, ttres[j].clone());
                    isout.addSlice(fpout);
                }

                ImagePlus ipout = new ImagePlus("tt3d,g="+IJ.d2s(gcsstd[i],1), isout);
                IJ.saveAsTiff(ipout, name);
                IJ.log("saving:\t"+name);

            }

        }
    }

    void init(float[] filtimg, int N, int M, int P, int[][] locxyz, float[] tness, int[] suppmap) {

        ArrayList[] Xall = new ArrayList[locxyz.length];

        for (int i = 0; i < Xall.length; i++) Xall[i] = new ArrayList<X>();

        int[] idx       = new int[locxyz.length];
        float[] score    = new float[locxyz.length];
        float[][] sumtness = new float[locxyz.length][1];

        for (int i = 0; i < locxyz.length; i++) {
            idx[i] = i;
            score[i] = corrAtLoc(filtimg, N, M, P, locxyz[i], tness, Xall[i], sumtness[i]);
            // highest corr Xk instance gets lowest index so that other directions are aligned with that direction
        }

        Tools.quicksort(score, idx); // set nr points per object here

        //*********************************************************
        ArrayList<X> Xsel           = new ArrayList<X>();
        ArrayList<Float> Xsel_CWS   = new ArrayList<Float>();

        X xj = new X(); // fun
        X x0 = new X(); // fun
        int nobjcount = 0;

        for (int i = locxyz.length-1; i >= 0; i--) { // go through the sorted values from the highest

            nobjcount++;

            if (nobjcount>nobjstart) break; // out of the loop

            for (int j = 0; j < Xall[idx[i]].size(); j++) {

                if (j==0) {
                    x0.set((X) Xall[idx[i]].get(j));
                    x0.w = x0.tness/sumtness[idx[i]][0];
                    Xsel.add(new X(x0));
                    if (Xsel_CWS.size()==0)     Xsel_CWS.add(x0.w);
                    else                        Xsel_CWS.add(x0.w + Xsel_CWS.get(Xsel_CWS.size()-1));
                }
                else {
                    xj.set((X) Xall[idx[i]].get(j));
                    xj.vx = x0.vx;
                    xj.vy = x0.vy;
                    xj.vz = x0.vz;
                    xj.w = xj.tness/sumtness[idx[i]][0];
                    Xsel.add(new X(xj));
                    Xsel_CWS.add(xj.w + Xsel_CWS.get(Xsel_CWS.size()-1));
                }

            }

        }

        //*********************************************************
        // resample with npcles
        npcles = nobjstart * ro;
        resample(Xsel, Xsel_CWS, npcles, Xk); // use selected guidepoints to sample initial PHD particles

        //*********************************************************
        // update supression map by setting step-sized neighbourhood to false
        for (int i = 0; i < Xk.size(); i++) {
            for (int j = 0; j < mm.psup.length; j++) {
                int x = Math.round(Xk.get(i).x) + mm.psup[j][0];
                int y = Math.round(Xk.get(i).y) + mm.psup[j][1];
                int z = Math.round(Xk.get(i).z) + mm.psup[j][2];
                if (x>=0 && x<N && y>=0 && y<M && z>=0 && z<P) {
                    if (suppmap[z*(N*M)+y*N+x]==-2) suppmap[z*(N*M)+y*N+x] = -1;
                }
            }
        }

        IJ.log("estimate init tness...");
        float[] tnessinitarray = new float[Xk.size()]; // store all tubularity values for initial particles

        for (int i = 0; i < Xk.size(); i++) tnessinitarray[i] = Xk.get(i).tness;

        tnessinit = Tools.median_Wirth(tnessinitarray);

        IJ.log("tnessinit="+IJ.d2s(tnessinit,2)+", clutter level="+IJ.d2s(cluttertness*tnessinit,2));

        String log = "init, ";
        log+="|X|="+Xk.size()+", ";
        phdmass = estimate(Xk, diam/2f, Yk);
        log+="PHD=" + IJ.d2s(phdmass,3)+", ";
        log+="|Y|=" + Yk.size() + ", ";
        npcles = Math.round(phdmass)*ro; // will be used to resample in next iteration
        log+="npcles=" + npcles+", ";
        Y.add((ArrayList<X>) Yk.clone());
        IJ.log(log);

    }

    static float gzx(X z, X x, float sigma) {
        return (float) Math.exp(-(Math.pow(x.x-z.x,2)+Math.pow(x.y-z.y,2)+Math.pow(x.z-z.z,2))/(2*Math.pow(sigma,2)));
    }

    private void locparticles(float locx, float locy, float locz, float dirx, float diry, float dirz, X[] atloc) {

        for (int i = 0; i < atloc.length; i++) {

            float x = (float) (locx+rndgen.nextGaussian()*0.5*step);
            float y = (float) (locy+rndgen.nextGaussian()*0.5*step);
            float z = (this.is2d)? 0 : (float) (locz+rndgen.nextGaussian()*0.5*step);

            atloc[i] = new X(x, y, z, dirx, diry, dirz, 1f, 1f, 1f/atloc.length, 1f); // Xk(x, y, z, vx, vy, vz, sig, corr, w, tness)

        }

    }

    private float corrAtLoc(float[] inimg, int N, int M, int P, int[] loc, float[] tness, ArrayList<X> atloc, float[] sumtness) {

        int[] dummy = new int[1];
        float score = 0;

        atloc.clear();
        sumtness[0] = 0;

        for (int i = 0; i < mm.cloudxyz.length; i++) { // it will be ro of them per location

            int x, y, z;
            float vx, vy, vz, corr, sig, tnesscurr;

            x = loc[0] + mm.cloudxyz[i][0];
            y = loc[1] + mm.cloudxyz[i][1];
            z = loc[2] + mm.cloudxyz[i][2];

            if (x>=0 && x<N && y>=0 && y<M && z>=0 && z<P) {

                vx = Float.NaN;
                vy = Float.NaN;
                vz = Float.NaN;

                corr = Float.NEGATIVE_INFINITY;
                sig = Float.NaN;

                for (int j = 0; j < vxyzUniform.length; j++) {

                    float curr_corr = zncc(x, y, z, vxyzUniform[j][0], vxyzUniform[j][1], vxyzUniform[j][2], inimg, N, M, P, dummy);
                    float curr_sig = gcsstd[dummy[0]];

                    if (curr_corr>corr) {
                        vx = vxyzUniform[j][0];
                        vy = vxyzUniform[j][1];
                        vz = vxyzUniform[j][2];
                        sig = curr_sig;
                        corr = curr_corr;
                    }

                }

                tnesscurr = tness[z*(N*M)+y*N+x];
                atloc.add(new X(x,y,z,      vx,vy,vz,       sig,corr,Float.NaN,tnesscurr));
                sumtness[0] += tnesscurr;
                score += corr;

            }
            else {
                return Float.NEGATIVE_INFINITY; // automatically discard the one that has pixels that are out
            }

        }

        X.sortByCorr(atloc);

        return score;

    }

    float zncc(float x, float y, float z, float vx, float vy, float vz, float[] img, int imgw, int imgh, int imgl, int[] gcsstdidx) {

        // loop through U,V,W offsets, sample from image by interpolating and store the values in img_vals[]
        // loop the same way as when templates were formed to be compliant
        float nrm = (float) Math.sqrt(vx*vx+vy*vy);

        if (nrm>0.001) {
            ux = ((vy<0)? -1 : 1)*(vy/nrm);
            uy = ((vy<0)?  1 :-1)*(vx/nrm);
            uz = 0;
        }
        else {
            ux = 1; uy = 0; uz = 0;
        }

        float ag = 0; // average

        if (this.is2d) { // one orthogonal is necessary

            wx = 0;
            wy = 0;
            wz = 0;

            for (int vv = -V2; vv <= V2; ++vv) { // template length: V*U*1
                for (int uu = -U2; uu <= U2; ++uu) {
                    //for (int ww = -W2; ww <= W2; ++ww) {

                    float xcoord = vv*(-vx) + uu*ux; // x components of all three orthogonals
                    float ycoord = vv*(-vy) + uu*uy;

                    if ( Math.floor(x+xcoord)<0   || Math.ceil(x+xcoord)>=imgw-1) return 0;
                    if ( Math.floor(y+ycoord)<0   || Math.ceil(y+ycoord)>=imgh-1) return 0;

                    float value = Tools.interp(x+xcoord, y+ycoord, 0, img, imgw, imgh, imgl);

                    img_vals[(vv+V2)*U+(uu+U2)] = value;
                    ag += value;

                    //}
                }
            }
        }
        else { // two orthogonals

            wx = uy*vz - uz*vy;
            wy = ux*vz - uz*vx;
            wz = ux*vy - uy*vx;

            for (int vv = -V2; vv <= V2; ++vv) { // template length: V*U*W
                for (int uu = -U2; uu <= U2; ++uu) {
                    for (int ww = -W2; ww <= W2; ++ww) {

                        float xcoord = vv*(-vx) + uu*ux + ww*wx;
                        float ycoord = vv*(-vy) + uu*uy + ww*wy;
                        float zcoord = vv*(-vz) + uu*uz + ww*wz;

                        if (Math.floor(x + xcoord)<0 || Math.ceil(x + xcoord)>=imgw-1) return 0;
                        if (Math.floor(y + ycoord)<0 || Math.ceil(y + ycoord)>=imgh-1) return 0;
                        if (Math.floor(z + zcoord)<0 || Math.ceil(z + zcoord)>=imgl-1) return 0;

                        float value = Tools.interp(x+xcoord, y+ycoord, z+zcoord, img, imgw, imgh, imgl);

                        img_vals[(vv+V2)*U*W+(ww+W2)*U+(uu+U2)] = value;
                        ag += value;

                    }
                }

            }

        }

        ag = ag / (U*V*W);

        // find correlation with corresponding template(s)
        float out_corr = Float.NEGATIVE_INFINITY;//-FLT_MAX; // ensure that at least one score will be higher
        int out_r_idx = -1; // index

        for (int tidx = 0; tidx < gcsstd_nr; ++tidx) {

            // check how it correlates with single[tidx], take max

            float corra = 0;
            float corrb = 0;
            float corrc = 0;

            for (int k = 0; k <U*W*V; k++) {
                corra += (img_vals[k]-ag) * (tt[tidx][k]-tta[tidx]);
                corrb += Math.pow(img_vals[k] - ag, 2);
                corrc += Math.pow(tt[tidx][k] - tta[tidx], 2);
            }

            float corr_val = (corrb*corrc>Float.MIN_VALUE)? (float) (corra / Math.sqrt(corrb * corrc)) : 0;

            if (corr_val>out_corr) {
                out_corr = corr_val;
                out_r_idx = tidx;// gcsstd[tidx];
            }

        }

        gcsstdidx[0] = out_r_idx;
        return out_corr; // compare out_corr later on

    }

    float zncc(X x, float[] img, int imgw, int imgh, int imgl, int[] gcsstdidx) {
        return zncc(x.x, x.y, x.z, x.vx, x.vy, x.vz, img, imgw, imgh, imgl, gcsstdidx);
    }

    float directionalavg(float[] locxyz, float[] dirxyz, float[] img, int iw, int ih, int il) {
        //
        float x  = locxyz[0];
        float y  = locxyz[1];
        float z  = locxyz[2];
        //
        float vx = dirxyz[0];
        float vy = dirxyz[1];
        float vz = dirxyz[2];

        // loop through U,V,W offsets, sample from image by interpolating and store the values in img_vals[]
        // loop the same way as when templates were formed to be compliant
        float nrm = (float) Math.sqrt(vx*vx+vy*vy);

        if (nrm>0.001) {
            ux = ((vy<0)? -1 : 1)*(vy/nrm);
            uy = ((vy<0)?  1 :-1)*(vx/nrm);
            uz = 0;
        }
        else {
            ux = 1; uy = 0; uz = 0;
        }

        float ag = 0; // average

        if (this.is2d) { // one orthogonal is necessary

            wx = 0;
            wy = 0;
            wz = 0;

            for (int vv = -V2; vv <= V2; ++vv) { // template length: V*U*1
                //for (int uu = -U2; uu <= U2; ++uu) {
                    //for (int ww = -W2; ww <= W2; ++ww) {

                    float uu = 0;
                    float xcoord = vv*(-vx) + uu*ux; // x components of all three orthogonals
                    float ycoord = vv*(-vy) + uu*uy;

                    if ( Math.floor(x+xcoord)<0   || Math.ceil(x+xcoord)>=iw-1) return 0;
                    if ( Math.floor(y+ycoord)<0   || Math.ceil(y+ycoord)>=ih-1) return 0;

                    float value = Tools.interp(x+xcoord, y+ycoord, 0, img, iw, ih, il);

//                    img_vals[(vv+V2)*U+(uu+U2)] = value;
                    ag += value;

                    //}
                //}
            }
        }
        else { // two orthogonals

            wx = uy*vz - uz*vy;
            wy = ux*vz - uz*vx;
            wz = ux*vy - uy*vx;

            for (int vv = -V2; vv <= V2; ++vv) { // template length: V*U*W
                //for (int uu = -U2; uu <= U2; ++uu) {
                    //for (int ww = -W2; ww <= W2; ++ww) {

                float uu = 0;
                float ww = 0;
                        float xcoord = vv*(-vx) + uu*ux + ww*wx;
                        float ycoord = vv*(-vy) + uu*uy + ww*wy;
                        float zcoord = vv*(-vz) + uu*uz + ww*wz;

                        if (Math.floor(x + xcoord)<0 || Math.ceil(x + xcoord)>=iw-1) return 0;
                        if (Math.floor(y + ycoord)<0 || Math.ceil(y + ycoord)>=ih-1) return 0;
                        if (Math.floor(z + zcoord)<0 || Math.ceil(z + zcoord)>=il-1) return 0;

                        float value = Tools.interp(x+xcoord, y+ycoord, z+zcoord, img, iw, ih, il);

//                        img_vals[(vv+V2)*U*W+(ww+W2)*U+(uu+U2)] = value;
                        ag += value;

                    //}
                //}

            }

        }

        return ag / (V);

    }

    boolean iter(int k, boolean initpred, float[] img, int N, int M, int P, float[] tness, int[] suppmap) {

        String log="k=" + k + ",\t|X|=" + Xk.size()+", PHD="+IJ.d2s(phdmass,2)+", ";

        //**************
        // prediction
        //**************
        int xi, yi, zi, xj, yj, zj, ii;
        double sc;
        float tnessval;
        X particle;

        XPk.clear();
        XPk_cws.clear();

        for (int i = 0; i < Xk.size(); i++) {

            X xp = Xk.get(i);

            xi = Math.round(xp.x);
            yi = Math.round(xp.y);
            zi = Math.round(xp.z);

            int ip = mm.getdirection(xp.vx, xp.vy, xp.vz); // pick one direction

            Arrays.fill(mm.pcws, 0); // cummulative weight sum for offsets to pick offset as random sample from
            // the distribution determined by p.d.f. that's multiplication of prior weights and the tness

            for (int j = 0; j < mm.p.length; j++) { // go through the Stepper offsets defined with step and kappa

                xj = xi + mm.p[j][0];
                yj = yi + mm.p[j][1];
                zj = zi + mm.p[j][2];
                ii = zj*(N*M)+yj*N+xj;

                if (xj>=0 && xj<N && yj>=0 && yj<M && zj>=0 && zj<P && suppmap[ii]==-2)
                    sc = mm.w[ip][j] * Math.pow(tness[ii],2);
                else
                    sc = 0;

                mm.pcws[j] = (float) ((j==0)? sc : (sc+mm.pcws[j-1]));

            }

            if (mm.pcws[mm.pcws.length-1]>Float.MIN_VALUE) { // sample from the cummulative distribution

                // *** sampling ni ***
                float wmass = mm.pcws[mm.pcws.length-1];
                float u1 = (wmass/ni)  * rndgen.nextFloat();

                int s = 0;

                for (int j = 0; j < ni; j++) {
                    float uj = u1 + j * (wmass/ni);
                    while (uj>mm.pcws[s]) s++;

                    // add s-th
                    xj = xi + mm.p[s][0];
                    yj = yi + mm.p[s][1];
                    zj = zi + mm.p[s][2];
                    ii = zj*(N*M)+yj*N+xj;
                    tnessval = tness[ii];

                    particle = new X(xj,yj,zj,   mm.u[s][0],mm.u[s][1],mm.u[s][2],    1,1, pS * xp.w * (1f/ni), tnessval);

                    XPk.add(particle);

                    if (XPk_cws.size()==0)  XPk_cws.add(particle.w);
                    else                    XPk_cws.add(particle.w + XPk_cws.get(XPk_cws.size()-1));

                }

            }

            if (initpred) {

                ip = mm.getdirection(-xp.vx, -xp.vy, -xp.vz);

                Arrays.fill(mm.pcws, 0);

                for (int j = 0; j < mm.p.length; j++) {

                    xj = xi + mm.p[j][0];
                    yj = yi + mm.p[j][1];
                    zj = zi + mm.p[j][2];
                    ii = zj*(N*M)+yj*N+xj;

                    if (xj>=0 && xj<N && yj>=0 && yj<M && zj>=0 && zj<P && suppmap[ii]==-2)
                        sc = mm.w[ip][j] * Math.pow(tness[ii],2);
                    else
                        sc = 0;

                    mm.pcws[j] = (float) ((j==0)? sc : (sc+mm.pcws[j-1]));

                }

                if (mm.pcws[mm.pcws.length-1]>Float.MIN_VALUE) { // sample from the cummulative distribution

                    // *** sampling ni ***
                    float wmass = mm.pcws[mm.pcws.length-1];
                    float u1 = (wmass/ni)  * rndgen.nextFloat();

                    int s = 0;

                    for (int j = 0; j < ni; j++) {
                        float uj = u1 + j * (wmass/ni);
                        while (uj>mm.pcws[s]) s++;

                        // add s-th
                        xj = xi + mm.p[s][0];
                        yj = yi + mm.p[s][1];
                        zj = zi + mm.p[s][2];
                        ii = zj*(N*M)+yj*N+xj;
                        tnessval = tness[ii];

                        particle = new X(xj,yj,zj,   mm.u[s][0],mm.u[s][1],mm.u[s][2],    1,1, pS * xp.w * (1f/ni), tnessval);

                        XPk.add(particle);

                        if (XPk_cws.size()==0)  XPk_cws.add(particle.w);
                        else                    XPk_cws.add(particle.w + XPk_cws.get(XPk_cws.size()-1));

                    }

                }
            } // initpred

        } // go through Xk PHD particles

        if (XPk.size()==0) {IJ.log("XPk.size()==0"); return false;}

        log+="|XP|=" + XPk.size()+", PHD="+IJ.d2s(XPk_cws.get(XPk_cws.size()-1),2)+", ";

        //**************
        // measure
        //**************
        ArrayList<Integer> cnt1 = new ArrayList<Integer>();
        ArrayList<X> Zk1 = new ArrayList<X>();
        measure(XPk, diam/2f, Zk1, cnt1); // be careful what to use here as dist

        // crop measurements so that they don't go into background
        Zk.clear();
        Zk_cnt.clear();

        // min values to be considered as a measurement
        int MIN_CNT = 2;
        float MIN_TNESS = (float) (0.01f*tnessinit);

        // select measurements - discard
        for (int i = 0; i < Zk1.size(); i++) {

            int x = Math.round(Zk1.get(i).x);
            int y = Math.round(Zk1.get(i).y);
            int z = Math.round(Zk1.get(i).z);
            Zk1.get(i).tness = tness[z*(N*M)+y*(N)+x]; // set tness value (it was NaN)

            if (Zk1.get(i).tness > MIN_TNESS && cnt1.get(i) >= MIN_CNT) { // it doesn't make sense to pick measurements below some level
                Zk.add(new X(Zk1.get(i)));
                Zk_cnt.add(cnt1.get(i));
            }

        }

        log+="|Z|="+Zk.size()+", ";

        if (Zk.size()==0) {IJ.log("Zk.size()==0"); return false;}

        //**************
        // update
        //**************
        XPk_cws = update(XPk, Zk, cnt1); // updated weights will be added to cws output

        phdmass = estimate(XPk, diam/2f, Yk); // final estimates

        if (Math.round(phdmass)<1) {IJ.log("Math.round(phdmass)<1"); return false;}

        if (Math.round(phdmass)>=OBJECT_LIMIT) {IJ.log("Math.round(phdmass)>=OBJECT_LIMIT, OBJECT_LIMIT=" + OBJECT_LIMIT); return false;}

        Y.add((ArrayList<X>) Yk.clone());
        log+="|Y|=" + Yk.size() + ", ";

        npcles = Math.round(phdmass)*ro;
        log+="npcles="+npcles;

        // resample
        resample(XPk, XPk_cws, npcles*((initpred)?2:1), Xk);

        //*********************************************************************
        // update supression map by setting step-sized neighbourhood to false
        for (int i = 0; i < Xk.size(); i++) {
            for (int j = 0; j < mm.psup.length; j++) {
                int x = Math.round(Xk.get(i).x) + mm.psup[j][0];
                int y = Math.round(Xk.get(i).y) + mm.psup[j][1];
                int z = Math.round(Xk.get(i).z) + mm.psup[j][2];
                if (x>=0 && x<N && y>=0 && y<M && z>=0 && z<P) {
                    if (suppmap[z*(N*M)+y*N+x]==-2) suppmap[z*(N*M)+y*N+x] = k;
                }
            }
        }

        IJ.log(log);

        return true;

    }

    private void measure(ArrayList<X> x, float dist, ArrayList<X> z, ArrayList<Integer> count) {

        float[] xw = new float[x.size()];

        for (int i = 0; i < x.size(); i++) {

            for (int j = 0; j < x.size(); j++) {
                float w =               x.get(j).w;
                float dx = x.get(i).x - x.get(j).x;
                float dy = x.get(i).y - x.get(j).y;
                float dz = x.get(i).z - x.get(j).z;
                float d2 = dx*dx + dy*dy + dz*dz;
                xw[i] += w * (float) Math.exp(-d2/(2*Math.pow(dist,2)));
            }

        }

        float[][] conv = Cluster.meanShift(x, xw, dist);
        int[] lab = Cluster.clustering(conv, dist);
        Cluster.extract(lab, conv, 1, Integer.MAX_VALUE, z, count); // extract() uses convergence values

    }

    private float estimate(ArrayList<X> x, float dist, ArrayList<X> y) {

        // cluster weighted phd particles using mean-shift and get the weighted mean out of each cluster
        float[] xw = new float[x.size()];
        float outcws = 0;
        for (int i = 0; i < x.size(); i++) {
            xw[i] = x.get(i).w;
            outcws += xw[i];
        }

        float[][] conv = Cluster.meanShift(x, xw, dist);
        int[] lab = Cluster.clustering(conv, dist);
        Cluster.estimate(lab, x, 1, Integer.MAX_VALUE, y); // estimate() uses x phd particles

        return outcws;
    }

    private ArrayList<Float> update(ArrayList<X> x, ArrayList<X> z, ArrayList<Integer> cnt) { // if updated with X instances, outputs CWS with updated weights, for resampling

        g = new float[x.size()][z.size()];
        Cz = new float[z.size()];

        // calculate Cz
        for (int j = 0; j < z.size(); j++) {

            // Cz[j] initialize with clutter PHD here
            Cz[j] = clutter(z.get(j).tness, kclutt, cluttertness*(float)tnessinit);

            for (int i = 0; i < x.size(); i++) {

                g[i][j] = gzx(z.get(j), x.get(i), (float)step);

                Cz[j] += pD * g[i][j] * x.get(i).w;

            }

        }
        
        // update weights using g, Cz and clutter phd
        ArrayList<Float> cws = new ArrayList<Float>(x.size());

        for (int i = 0; i < x.size(); i++) {

            float wup = (1-pD) * x.get(i).w; // wup == updated weight for i-th particle

            for (int j = 0; j < z.size(); j++) {
                wup += (pD*g[i][j]*x.get(i).w)/Cz[j]; // update particle weight equation
            }

            x.get(i).w = wup;

            if (i==0)   cws.add(i, wup);
            else        cws.add(i, wup + cws.get(cws.size()-1));

        }

        return cws;

    }

    public static int[] cluster(ArrayList<X> Xlist) {

        int[] labels = new int[Xlist.size()];
        for (int i = 0; i < labels.length; i++) labels[i] = i;

        for (int i = 0; i < Xlist.size(); i++) {

            // one versus the rest
            for (int j = 0; j < Xlist.size(); j++) {

                if (i!=j) {

                    double dst2 	=
                                    Math.pow(Xlist.get(i).x-Xlist.get(j).x, 2) +
                                    Math.pow(Xlist.get(i).y-Xlist.get(j).y, 2) +
                                    Math.pow(Xlist.get(i).z-Xlist.get(j).z, 2);

                    double rd2 		=
                                    Math.pow(Xlist.get(i).sig+Xlist.get(j).sig, 2);


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

    public static ArrayList<Z> extract(int[] labels, ArrayList<X> vals) { //int[] vals

        boolean[] checked = new boolean[labels.length];
        ArrayList<Z> out = new ArrayList<Z>();

        for (int i = 0; i < labels.length; i++) {
            if (!checked[i]) {

                X cc = vals.get(i); //[ i ]; // idxs[i]
                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < labels.length; j++) {
                    if (!checked[j]) {
                        if (labels[j]==labels[i]) {

                            cc.x += vals.get(j).x;
                            cc.y += vals.get(j).y;
                            cc.z += vals.get(j).z;
                            cc.sig += vals.get(j).sig;

                            count++;
                            checked[j] = true;

                        }
                    }
                }

                out.add(new Z(cc.x/count, cc.y/count, cc.z/count, cc.sig/count, count));
//                out.add(new float[]{centroid/count, count});

            }
        }

        return out;

    }

    private ArrayList<X> importsamp(ArrayList<X> xlist, float[] cws, int N) {
        // systematic resampling, Beyond Kalman Filtering, Ristic et al.
        // xlist elements will be randomly sampled using the weights from cws[]

        float totalmass = cws[cws.length-1]; // total cummulative mass

        float u1 = (totalmass/(float)N) * rndgen.nextFloat();

        ArrayList<X> out = new ArrayList<X>(N);
        out.clear();
        int i = 0;
        for (int j = 0; j < N; j++) {
            float uj = u1 + j*(totalmass/(float)N);
            while (uj > cws[i]) i++;
            X tt = new X(xlist.get(i));
            tt.w = totalmass/(float)N;
            out.add(tt);
        }

        return out;

    }

    void resample(ArrayList<X> Xin, ArrayList<Float> Cws, int N, ArrayList<X> Xout) {
        // systematic resampling algorithm, Beyond Kalman Filtering, Ristic et al.
        float wmass = Cws.get(Cws.size()-1);
        float u1 = (wmass/N) * rndgen.nextFloat();

        Xout.clear();
        
        int i = 0;

        for (int j = 0; j < N; j++) {
            float uj = u1 + j*(wmass/(float)N);
            while (uj > Cws.get(i)) i++;
            X tt = new X(Xin.get(i));
            tt.w = wmass/N;
            Xout.add(tt);
        }

    }

//    public static float clutter(int c, float K, int cc) {
//        return (float) Math.exp(-K*(c-cc));
//    }

    public static float clutter(float t, float K, float tc) {
        return (float) Math.exp(-K*(t-tc));
    }

}
