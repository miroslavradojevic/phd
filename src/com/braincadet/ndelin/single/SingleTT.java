package com.braincadet.ndelin.single;

/*
Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
Created by miroslav on 9/21/15.
*/

import com.braincadet.ndelin.fun.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.Arrays;
import java.util.Random;

public class SingleTT {

    Random rndgen = new Random();

    static float gcsstd_min  = 1f;
    static float gcsstd2rad  = 2.5f;
    static float gcsstd_step = 1f;
    static float step = 2f;  // object tracking step
    static float locstddev = 0.25f;
    static int   ndirs2d = 15; // 1/2 circle, radius=step
    static int   ndirs3d = 25; // 1/2 sphere, radius=step
    static float ratioNthr = 0.5f;
    static float K = 4; // used in likelihood calculation exp(K*zncc)

    boolean   is2d;
    float     gcsstd_max;
    int       gcsstd_nr;
    float[]   gcsstd;
    int       sc;
    int       U, W, V, U2, W2, V2;
    float[]   img_vals;
    float[][] tt;
    float[]   tta;

    // for orthogonals in znccX
    float       ux, uy, uz;
    float       wx, wy, wz;

    //used in geometry transformations, affine
    float[][]   R33     = new float[3][3];
    float[][]   R22     = new float[2][2];

    float[]     v13    = new float[3];
    float[]     v12    = new float[2];

    public SingleTT(int sc, boolean is2d) {

        this.sc = sc;
        this.is2d = is2d;

        // ------------------------------------------------------------------------
        // templates formation (templates used to calculate the likelihood)
        // template size U*W*V, W is added for 3d (it is 1 in 2d), V is longitudinal
        // ------------------------------------------------------------------------

        U2 = this.sc/2;
        U = 2*U2 + 1; // 1st orthogonal

        if (this.is2d) {
            W2 = 0;
            W = 1;
        }
        else {
            W2 = this.sc/2;
            W = 2*W2 + 1;
        }

        V2 = this.sc/2;//(this.is2d)?sc/2:1;
        V = 2*V2 + 1;

        // gcsstd define
        gcsstd_max = (U2/gcsstd2rad);
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

    void demo(int np, int ni, float angrad, String outputdir) {

        PrintWriter logWriter = null;

        if (outputdir!=null) {

            String outfile = outputdir + File.separator + "demo,sc=" + IJ.d2s(sc,1) + ",np=" + IJ.d2s(np,0) + ",ni=" + IJ.d2s(ni,0) + ".swc";

            try {
                logWriter = new PrintWriter(outfile);
                logWriter.print("");
                logWriter.close();
            } catch (FileNotFoundException ex) {
            }

            try {
                logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));
                logWriter.println("#demo()");
            } catch (IOException e) {
            }
        }

        float vx    = 1f;
        float vy    = 0f;
        float vz    = 0f;

        float[][][] pxyz = new float[ni][np][3];
        float[][][] vxyz = new float[ni][np][3];
        float[][]   uxyz;

        logWriter.println("##### 2d #####");

        int swci = 0;
        logWriter.println((++swci) + " " + 1 + " " + 0 + " " + 0 + " " + 0 + " " + 0.15 + " " + (-1));
        int swcroot= swci;

        uxyz = Tools.calcdirs180(true, ndirs2d); // unit directions half circle (0,0,0)
        for (int j = 0; j < uxyz.length; j++) {
            logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(uxyz[j][0], 3) + " " + IJ.d2s(uxyz[j][1], 3) + " " + IJ.d2s(uxyz[j][2], 3) + " " + 0.1 + " " + swcroot);
            logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(-uxyz[j][0], 3) + " " + IJ.d2s(-uxyz[j][1], 3) + " " + IJ.d2s(-uxyz[j][2], 3) + " " + 0.1 + " " + swcroot);
        }

        logWriter.println("########");

        logWriter.println((++swci) + " " + 1 + " " + 10 + " " + 0 + " " + 0 + " " + 0.1 + " " + (-1));
        swcroot = swci;
        logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(10 + sc * vx, 3) + " " + IJ.d2s(0 + sc * vy, 3) + " " + IJ.d2s(0 + sc * vz, 3) + " " + 0.1 + " " + swcroot);

        // _init iteration
        for (int i = 0; i < np; i++) { // random initial locations, start with (10,0,0) for visualization
            pxyz[0][i][0] = (float) (10+rndgen.nextGaussian()*locstddev);
            pxyz[0][i][1] = (float) (0 +rndgen.nextGaussian()*locstddev);
            pxyz[0][i][2] = 0;
            vxyz[0][i][0] = vx;
            vxyz[0][i][1] = vy;
            vxyz[0][i][2] = vz;
        }

        boolean is2dprev = this.is2d;
        this.is2d = true; // for followup()
        for (int i = 1; i < ni; i++) {
            for (int k = 0; k < np; k++) {
                followup(pxyz[i-1][k][0],pxyz[i-1][k][1], pxyz[i-1][k][2], vxyz[i-1][k][0], vxyz[i-1][k][1], vxyz[i-1][k][2], angrad, pxyz[i][k], vxyz[i][k]);
            }
        }

        for (int i = 0; i < ni; i++) {
            for (int k = 0; k < np; k++) {
                logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(pxyz[i][k][0], 3) + " " + IJ.d2s(pxyz[i][k][1], 3) + " " + IJ.d2s(pxyz[i][k][2], 3) + " " + 0.1 + " " + (-1));
            }
        }

        logWriter.println("##### 3d #####");

        logWriter.println((++swci) + " " + 1 + " " + 0 + " " + 30 + " " + 0 + " " + 0.1 + " " + (-1));
        swcroot= swci;

        uxyz = Tools.calcdirs180(false, ndirs3d); // unit directions half sphere
        for (int j = 0; j < uxyz.length; j++) {
            logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(0 + uxyz[j][0], 3) + " " + IJ.d2s(30 + uxyz[j][1], 3) + " " + IJ.d2s(0 + uxyz[j][2], 3) + " " + 0.1 + " " + swcroot);
            logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(0 - uxyz[j][0], 3) + " " + IJ.d2s(30 - uxyz[j][1], 3) + " " + IJ.d2s(0 - uxyz[j][2], 3) + " " + 0.1 + " " + swcroot);
        }

        logWriter.println((++swci) + " " + 1 + " " + 10 + " " + 30 + " " + 0 + " " + 0.1 + " " + (-1));
        swcroot = swci;
        logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(10 + sc * vx, 3) + " " + IJ.d2s(30 + sc * vy, 3) + " " + IJ.d2s(0 + sc * vz, 3) + " " + 0.1 + " " + swcroot);

        for (int i = 0; i < np; i++) { // random initial locations, start with (10,0,0)
            pxyz[0][i][0] = (float) (10+rndgen.nextGaussian()*locstddev);
            pxyz[0][i][1] = (float) (30+rndgen.nextGaussian()*locstddev);
            pxyz[0][i][2] = (float) (0 +rndgen.nextGaussian()*locstddev);
            vxyz[0][i][0] = vx;
            vxyz[0][i][1] = vy;
            vxyz[0][i][2] = vz;
        }

        this.is2d = false;
        for (int i = 1; i < ni; i++) {
            for (int k = 0; k < np; k++) {
                followup(pxyz[i-1][k][0],pxyz[i-1][k][1], pxyz[i-1][k][2], vxyz[i-1][k][0], vxyz[i-1][k][1], vxyz[i-1][k][2], angrad, pxyz[i][k], vxyz[i][k]);
            }
        }

        for (int i = 0; i < ni; i++) {
            for (int k = 0; k < np; k++) {
                logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(pxyz[i][k][0], 3) + " " + IJ.d2s(pxyz[i][k][1], 3) + " " + IJ.d2s(pxyz[i][k][2], 3) + " " + 0.1 + " " + (-1));
            }
        }

        logWriter.close();
        this.is2d = is2dprev;

    }

    float[] locdirection(
            float[] inimg,
            int img_w,
            int img_h,
            int img_l,
            float atX,
            float atY,
            float atZ,
            String outputdir) {

        float[][] vxyz = Tools.calcdirs180(this.is2d, (this.is2d ? ndirs2d : ndirs3d)); // store directions

        // locations (20) with gaussian random with locsstddev
        float[][] locsxyz = new float[20][3];
        for (int i = 0; i < locsxyz.length; i++) { // random initial locations, start with (10,0,0)
            locsxyz[i][0] = (float) (atX+rndgen.nextGaussian()*locstddev);
            locsxyz[i][1] = (float) (atY+rndgen.nextGaussian()*locstddev);
            locsxyz[i][2] = (this.is2d)? 0 : (float) (atZ+rndgen.nextGaussian()*locstddev);
        }

        // loop through the set of directions and take the one with the highest average correlation
        float maxcorr = Float.NEGATIVE_INFINITY;
        int diridx = -1;
        int[] dummy = new int[1];
        for (int i = 0; i < vxyz.length; i++) {

            float avgcorr = 0;

            for (int j = 0; j < locsxyz.length; j++) {
                avgcorr += znccX(locsxyz[j][0], locsxyz[j][1], locsxyz[j][2], vxyz[i][0], vxyz[i][1], vxyz[i][2], inimg, img_w, img_h, img_l, dummy);
            }

            avgcorr /= (float) locsxyz.length;

            if (avgcorr>maxcorr) {
                maxcorr = avgcorr;
                diridx = i;
            }

        }

        IJ.log("at " + Arrays.toString(vxyz[diridx]) + " -> " + IJ.d2s(maxcorr,2));

        if (outputdir!=null) {

            String outfile = outputdir+File.separator+"locdirection,x="+IJ.d2s(atX,3)+",y="+IJ.d2s(atY,3)+",z="+IJ.d2s(atZ,3)+".swc";

            PrintWriter logWriter = null;
            try {
                logWriter = new PrintWriter(outfile);
                logWriter.print("");
                logWriter.close();
            } catch (FileNotFoundException ex) {}

            try {
                logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));
                logWriter.println("#locdirection()");
            } catch (IOException e) {}

            int swci = 0;
            for (int i = 0; i < locsxyz.length; i++) {
                if (logWriter != null) {
                    logWriter.println((++swci) + " " + 2 + " " + IJ.d2s(locsxyz[i][0],3) + " " + IJ.d2s(locsxyz[i][1],3) + " " + IJ.d2s(locsxyz[i][2],3) + " " + 1 + " " + -1);
                }
            }
            if (logWriter != null) {
                logWriter.println((++swci) + " " + 3 + " " + IJ.d2s(atX,3) + " " + IJ.d2s(atY,3) + " " + IJ.d2s(atZ,3) + " " + 1 + " " + -1);
            }
            int swcroot = swci;
            if (logWriter != null) {
                logWriter.println((++swci) + " " + 3 + " " + IJ.d2s(atX + sc * vxyz[diridx][0],3) + " " + IJ.d2s(atY + sc * vxyz[diridx][1],3) + " " + IJ.d2s(atZ + sc * vxyz[diridx][2],3) + " " + 1 + " " + swcroot);
            }
            logWriter.close();
            IJ.log("exported\t"+outfile);

        }

        return vxyz[diridx].clone();

    }

    void track(
            float[] inimg,
            int img_w,
            int img_h,
            int img_l,
            float atX,
            float atY,
            float atZ,
            float atVX,
            float atVY,
            float atVZ,
            int sc, float angrad, int ni, int np, String outputdir) {

        // track logger initialize
        PrintWriter logWriter = null;
        String outfile = outputdir + File.separator + "track,x=" + IJ.d2s(atX,2) + ",y=" + IJ.d2s(atY,2) + ",z=" + IJ.d2s(atZ,2) + ".swc";
        IJ.log(outfile);
        if (outputdir != null) {

            try {
                logWriter = new PrintWriter(outfile);
                logWriter.print("#track()");
                logWriter.close();
            } catch (FileNotFoundException ex) {
            }

            try {
                logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));
                logWriter.println("");
            } catch (IOException e) {
            }

            logWriter.println("#(x,y,z)=" + atX + ", " + atY + ", " + atZ);
            logWriter.println("#(vx,vy,vz)=" + atVX + ", " + atVY + ", " + atVZ);
            logWriter.println("#sc=" + sc);
            logWriter.println("#angrad=" + angrad);
            logWriter.println("#ni=" + ni);
            logWriter.println("#np=" + np);

        }
        /////////////////////////////////////////////////////////////////////////////////////
        float[][] x = new float[np][3];     // particles x_k^i
        float[][] vx = new float[np][3];    // particle directions

        float[] xest = new float[3];        // state estimation
        float[] w = new float[np];          // particle weights w_k^i
        float[] wprev = new float[np];      // weights from the previous iteration
        float[] g = new float[np];          // lhoods at iteration
        int[] gcsstdidx = new int[1];       // index of the gaussian cross section std.
        int Nthr = (int) Math.floor(ratioNthr*np); Nthr = (Nthr<1)? 1 : Nthr;
        float[] c = new float[np];          // csw, for resampling

        float[][] xr = new float[np][3];    // to store resampled states (locations and directions)
        float[][] vxr = new float[np][3];

        /////////////////////////////////////////////////////////////////////////////////////
        // smc-pf implementation
        // initial particles, weights and estimate
        Arrays.fill(xest, 0f);
        int iter = 0;
        for (int i = 0; i < np; i++) {

            x[i][0] = (float) (atX + rndgen.nextGaussian()*locstddev);
            x[i][1] = (float) (atY + rndgen.nextGaussian()*locstddev);
            x[i][2] = (float) (atZ + rndgen.nextGaussian()*locstddev);

            vx[i][0] = atVX;
            vx[i][1] = atVY;
            vx[i][2] = atVZ;

            w[i] = 1f/(float)np;

            xest[0] += w[i] * x[i][0];
            xest[1] += w[i] * x[i][1];
            xest[2] += w[i] * x[i][2];

        }

        iter++;

        if (outputdir != null) {
            logWriter.println(iter + " " + 3 + " " + IJ.d2s(xest[0],3) + " " + IJ.d2s(xest[1],3) + " " + IJ.d2s(xest[2],3) + " " + 1 + " " + -1);
        }

        while (iter<ni) {

            float w1 = 0;

            for (int i = 0; i < np; i++) {

                float xprev = x[i][0];
                float yprev = x[i][1];
                float zprev = x[i][2];

                float vxprev = vx[i][0];
                float vyprev = vx[i][1];
                float vzprev = vx[i][2];

                wprev[i] = w[i];

                w[i] = followup(
                        xprev,   yprev,   zprev,
                        vxprev,  vyprev,  vzprev,
                        angrad,
                        x[i],
                        vx[i]); // prediction

                w1 += w[i];

                g[i] = znccX(
                        x[i][0],     x[i][1],     x[i][2],
                        vx[i][0],    vx[i][1],    vx[i][2],
                        inimg, img_w, img_h, img_l, gcsstdidx); // likelihood

                g[i] = (float) Math.exp(K*g[i]); // g[i] in [-1, 1]

            }

            float w2 = 0;

            // update - compute weights
            for (int i = 0; i < np; i++) {
                w[i] = wprev[i] * g[i] * (w[i]/w1);
                w2 += w[i];
            }

            // normalize weights & calculate estimate before resampling & calculate effective sample size
            Arrays.fill(xest, 0f);
            float Neff = 0;
            for (int i = 0; i < np; i++) {
                w[i] /= w2;
                xest[0] += w[i] * x[i][0];
                xest[1] += w[i] * x[i][1];
                xest[2] += w[i] * x[i][2];
                Neff += w[i] * w[i];

            }

            Neff = 1f/Neff;

            IJ.log("iter0 = " + iter);
//            IJ.log("iter0=" + iter0 + ", Neff=" + Neff + ", g25= " + quantile(g, 25, 100) + ", g50=" + quantile(g, 50, 100) + ", g75=" + quantile(g, 75, 100) + ", th=" + Math.exp(K*0.5f));
//            IJ.log("---");
//            IJ.log(Arrays.toString(g));
//            IJ.log("---");

            if (outputdir != null) {
                logWriter.println(iter + " " + 3 + " " + IJ.d2s(xest[0],3) + " " + IJ.d2s(xest[1],3) + " " + IJ.d2s(xest[2],3) + " " + 1 + " " + (iter-1));
            }

            if (Neff<Nthr) {

                // systematic resampling algorithm, from Beyond Kalman Filtering, Ristic et al.
                c[0] = w[0];
                for (int i = 1; i < np; i++) {
                    c[i] = c[i-1] + w[i];
                }

                int i = 0;

                float u1 = (1f/(float)np) * rndgen.nextFloat(); // nextFloat() u[0,1) -> u[0, 1/np), it was supposed to be [0, 1/np]

                for (int j = 0; j < np; j++) {

                    float uj = u1 + j/(float)np;

                    while (uj > c[i]) {
                        i++;
                    }

                    xr[j][0] = (float) (x[i][0] + rndgen.nextGaussian()*locstddev);
                    xr[j][1] = (float) (x[i][1] + rndgen.nextGaussian()*locstddev);
                    xr[j][2] = (float) (x[i][2] + rndgen.nextGaussian()*locstddev);
                    vxr[j][0] = vx[i][0];
                    vxr[j][1] = vx[i][1];
                    vxr[j][2] = vx[i][2];

                }

                // assign resampled values to x states and reset the weights
                for (int k = 0; k < np; k++) {
                    x[k][0] = xr[k][0];
                    x[k][1] = xr[k][1];
                    x[k][2] = xr[k][2];
                    vx[k][0] = vxr[k][0];
                    vx[k][1] = vxr[k][1];
                    vx[k][2] = vxr[k][2];
                    w[k] = 1f/(float)np;
                }

            }

            iter++;

        }

        /////////////////////////////////////////////////////////////////////////////////////
        if (outputdir != null) {
            logWriter.close();
            IJ.log("exported\t"+outfile);
        }

    }

    void exporttemplates(String outdir) {

        for (int i = 0; i < tt.length; i++) {

            String name = outdir+ File.separator+"template,g="+IJ.d2s(gcsstd[i],1)+".tif";

            if (this.is2d) {
                ImageStack isout = new ImageStack(U,V); // W=1
                FloatProcessor fpout = new FloatProcessor(U, V, tt[i].clone());
                isout.addSlice(fpout);
                ImagePlus ipout = new ImagePlus("g="+IJ.d2s(gcsstd[i],1), isout);
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

                ImagePlus ipout = new ImagePlus("g="+IJ.d2s(gcsstd[i],1), isout);
                IJ.saveAsTiff(ipout, name);
                IJ.log("saving:\t"+name);

            }

        }
    }

    float znccX(float x, float y, float z,
                float vx, float vy, float vz,
                float[] img, int img_w, int img_h, int img_l,
                int[] gcsstd_idx  // side-output
    ){

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

                        if ( Math.floor(x+xcoord)<0   || Math.ceil(x+xcoord)>=img_w-1) return 0;
                        if ( Math.floor(y+ycoord)<0   || Math.ceil(y+ycoord)>=img_h-1) return 0;

                        float value = Tools.interp(x+xcoord, y+ycoord, 0, img, img_w, img_h, img_l);

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

                        if (Math.floor(x + xcoord)<0 || Math.ceil(x + xcoord)>=img_w-1) return 0;
                        if (Math.floor(y + ycoord)<0 || Math.ceil(y + ycoord)>=img_h-1) return 0;
                        if (Math.floor(z + zcoord)<0 || Math.ceil(z + zcoord)>=img_l-1) return 0;

                        float value = Tools.interp(x+xcoord, y+ycoord, z+zcoord, img, img_w, img_h, img_l);

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

            float corr_val = (corrb*corrc>0.0001)? (float) (corra / Math.sqrt(corrb * corrc)) : 0;

            if (corr_val>out_corr) {
                out_corr = corr_val;
                out_r_idx = tidx;// gcsstd[tidx];
            }

        }

        gcsstd_idx[0] = out_r_idx;
        return out_corr;

    }

    float followup(float x, float y, float z, float vx, float vy, float vz, float alfarad, float[] x1, float[] vx1) {

        if (this.is2d) { // z=0, vz=0

            float theta = (float) (rndgen.nextGaussian()*alfarad); // gaussian

            float ux = (float) Math.sin(theta);
            float uy = (float) Math.cos(theta);

            Tools.rotation_matrix(0, 1, vx, vy, R22);
            Tools.rotation_apply(R22, ux, uy, v12);

            x1[0] = x + v12[0] * this.step;
            x1[1] = y + v12[1] * this.step;
            x1[2] = 0;

            vx1[0] = v12[0];
            vx1[1] = v12[1];
            vx1[2] = 0;

            return uy;

        }
        else { // z!=0, vz!=0

            float phi = rndgen.nextFloat()*2f*3.14f; // uniform
            float theta = (float) (rndgen.nextGaussian()*alfarad); // gaussian

            float ux = (float) (Math.sin(theta) * Math.cos(phi));
            float uy = (float) (Math.sin(theta) * Math.sin(phi));
            float uz = (float)  Math.cos(theta);

            Tools.rotation_matrix(0, 0, 1, vx, vy, vz, R33); // R33 aligns (0,0,1) to (vx,vy,vz)
            Tools.rotation_apply(R33, ux, uy, uz, v13); // rotate using R33, store output in v13[]

            x1[0] = x + v13[0] * this.step;
            x1[1] = y + v13[1] * this.step;
            x1[2] = z + v13[2] * this.step; //((this.is2d)? 0 : 1);

            vx1[0] = v13[0];
            vx1[1] = v13[1];
            vx1[2] = v13[2];

            return uz;

        }

    }

}
