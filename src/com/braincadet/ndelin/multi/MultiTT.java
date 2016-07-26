package com.braincadet.ndelin.multi;

import com.braincadet.ndelin.fun.Tools;
import com.braincadet.ndelin.swc.Node;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.*;

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

    int         diam = 10;                  // tube diameter in pixels, fixed value, used only for the initial direction
    float       kernel_radius = -1;
    int         step = -1;                  // radius for sampling new particles

    float       kappa = Float.NaN;          // von Mises distribution kappa = 0.5, 1, 2
    float       pS = Float.NaN;             // probability
    float       pD = Float.NaN;             // detection probability
    public float cluttertness = Float.NaN;  // 0-background, 1-tubularity level at which _init objects are
    float       kclutt = Float.NaN;         // drop in PHD measure of the clutter
//    double      tnessinit = Double.NaN;   // tubularity level at initial objects

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
//    public ArrayList<ArrayList<Integer>> XPCk;  // XPk list indexes per cluster (to select those that will be used in Y and suppmap)
    public ArrayList<Float> XPk_cws;            // cummulative weight sum for the predicted particles

    public ArrayList<X> Zk;                     // measurements
    public ArrayList<X> ZPk;                    // predicted paritcles used to get the measurements
    public ArrayList<ArrayList<Integer>> ZPCk;  // ZPk list indexes per cluster (to select those that will be used in Zk)

    public ArrayList<Node> Y;

    public static int   OBJECT_LIMIT = 200;     // limit amount of objects filtered
    public static int   MAXITER     = Integer.MAX_VALUE;
    public static float EPSILON2    = 0.000001f;

//    private X[]         Xpred;   // auxilliary storage for predicted particles
    public int          npcles;  // number of particles approximating the probability density
    public float        phdmass; //

    public float[][]    g;
    public float[]      Cz;

    public static int MAX_RADIUS = 15;
    public static int[][][] offxyz;

    float[]         xw;
    float[][]       conv;
    int[]           labels;
    boolean[]       checked;
    ArrayList[]     nbridxs;

    public boolean verbose = false; // to show the particle log
    public int R_supp = 0;
    public float gzx_sigma = 2f; // used in calculating likelihood, the distance towards measurement
    public static float weight_deg = 5;
    public static float weight_deg77 = 2; // for init locations
    int MIN_CLUST_SIZE = -1; //(int) Math.round(0.1*x.size()); will refer to ni
    float wmin = 0.4f;

    public MultiTT(boolean is2d, int no, int ro, int ni, int krad, int step, float kappa, float pS, float pD, float cluttertness, float kclutt) {

        offxyz = new int[MAX_RADIUS+1][][];
        offxyz[0] = new int[1][3];
        offxyz[0][0][0] = 0;
        offxyz[0][0][1] = 0;
        offxyz[0][0][2] = 0;

        for (int r = 1; r <= MAX_RADIUS; r++) {

            ArrayList<Integer> px = new ArrayList<Integer>();
            ArrayList<Integer> py = new ArrayList<Integer>();
            ArrayList<Integer> pz = new ArrayList<Integer>();

            for (int x = -r; x <= r; x++) {
                for (int y = -r; y <= r; y++) {
                    if (is2d) {
                        if (x*x+y*y<=r*r) {
                            px.add(x); py.add(y); pz.add(0);
                        }
                    }
                    else {
                        for (int z = -r; z <= r; z++) {
                            if (x*x+y*y+z*z<=r*r) {
                                px.add(x); py.add(y); pz.add(z);
                            }
                        }
                    }
                }
            }

//            IJ.log(px.size() + " points in the neighbourhood for R="+r);

            offxyz[r] = new int[px.size()][3];

            for (int i = 0; i < px.size(); i++) {
                offxyz[r][i][0] = px.get(i);
                offxyz[r][i][1] = py.get(i);
                offxyz[r][i][2] = pz.get(i);
            }

        }

        this.is2d = is2d;
        this.nobjstart = (no>OBJECT_LIMIT)?OBJECT_LIMIT:no;
        this.ro = ro;
        this.ni = ni;
        this.MIN_CLUST_SIZE = 2*ni;
        this.kernel_radius = krad;
        this.step = step;
        this.kappa = kappa;
        this.pS = pS;
        this.pD = pD;
        this.cluttertness = cluttertness;
        this.kclutt = kclutt;

        mm = new Stepper(this.step, this.is2d, this.kappa, this.ro);
//        IJ.log("TOTAL MEMORY = "+ MTracker.bytesToMegabytes( Runtime.getRuntime().totalMemory()) + " Mb" );

        Xk = new ArrayList<X>();
        XPk = new ArrayList<X>();
//        XPCk = new ArrayList<ArrayList<Integer>>();
        XPk_cws = new ArrayList<Float>();

        Zk = new ArrayList<X>();
        ZPk = new ArrayList<X>();
        ZPCk = new ArrayList<ArrayList<Integer>>();

        Y = new ArrayList<Node>();
        Y.add(null);

        // variables used in calculating iterations
        int CAPACITY = OBJECT_LIMIT*ro*ni;
        xw      = new float[CAPACITY];
        conv    = new float[CAPACITY][3];
        labels  = new int[CAPACITY];
        checked = new boolean[CAPACITY];
        nbridxs = new ArrayList[CAPACITY];
        for (int i = 0; i < CAPACITY; i++)
            nbridxs[i] = new ArrayList<Integer>();

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
        V2 = this.diam/2;
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

    public void exporttemplates(String outdir) {

        for (int i = 0; i < tt.length; i++) {

            String name = outdir+ File.separator+"template,g="+IJ.d2s(gcsstd[i],1)+".tif";

            if (this.is2d) {
                ImageStack isout = new ImageStack(U,V); // W=1
                FloatProcessor fpout = new FloatProcessor(U, V, tt[i].clone());
                isout.addSlice(fpout);
                ImagePlus ipout = new ImagePlus("tt2d,g="+IJ.d2s(gcsstd[i],1), isout);
                IJ.saveAsTiff(ipout, name);

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

            }

        }
    }

    static float gzx(X z, X x, float sigma) {
        return (float) Math.exp(-(Math.pow(x.x-z.x,2)+Math.pow(x.y-z.y,2)+Math.pow(x.z-z.z,2))/(2*Math.pow(sigma,2)));
    }

    private void locparticles(float locx, float locy, float locz, float dirx, float diry, float dirz, X[] atloc) {

        for (int i = 0; i < atloc.length; i++) {

            float x = (float) (locx+rndgen.nextGaussian()*0.5*step);
            float y = (float) (locy+rndgen.nextGaussian()*0.5*step);
            float z = (this.is2d)? 0 : (float) (locz+rndgen.nextGaussian()*0.5*step);

            atloc[i] = new X(x, y, z, dirx, diry, dirz, 1f/atloc.length, 1f); // Xk(x, y, z, vx, vy, vz, sig, corr, w, tness) 1f, 1f,

        }

    }

    private boolean extractloc(float[] inimg, int N, int M, int P, int[] loc, float[] tness, ArrayList<X> atloc) {

        int[] dummy = new int[1];
        atloc.clear();
        float sumw = 0;

        for (int i = 0; i < mm.cloudxyz.length; i++) { // ro particles per location

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
//                sig = Float.NaN;

                for (int j = 0; j < vxyzUniform.length; j++) {

                    float curr_corr = zncc(x,y,z, vxyzUniform[j][0], vxyzUniform[j][1], vxyzUniform[j][2], inimg, N, M, P, dummy);
                    float curr_sig = gcsstd[dummy[0]];

                    if (curr_corr>corr) {
                        vx = vxyzUniform[j][0];
                        vy = vxyzUniform[j][1];
                        vz = vxyzUniform[j][2];
//                        sig = curr_sig;
                        corr = curr_corr;
                    }

                }

                tnesscurr = tness[z*(N*M)+y*N+x];
                float w = tnesscurr; // tnesscurr  corr
                atloc.add(new X(x,y,z,      vx,vy,vz,        w,tnesscurr));//sig,corr,
                sumw += w;

            }
            else {
                return false;
            }

        }

//        X.sortByCorr(atloc); // not necessary any more

        for (int i = 0; i < atloc.size(); i++) {

            if (sumw>Float.MIN_VALUE) {
                atloc.get(i).w /= sumw;
            }
            else {
                atloc.get(i).w = 1f/atloc.size();
            }

        }

        return true;

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
                atloc.add(new X(x,y,z,      vx,vy,vz,       Float.NaN,tnesscurr)); // sig,corr,
                sumtness[0] += tnesscurr;
                score += corr;

            }
            else {
                return Float.NEGATIVE_INFINITY; // automatically discard the one that has pixels that are out
            }

        }

//        X.sortByCorr(atloc); // not necessary

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

    public void _init(int n, int radius, ArrayList<Integer> l, ArrayList<Float> w, ArrayList<int[]> N_o, float[] img, int N, int M, int P, float[] tness, int[] suppmap) {

        ArrayList<Integer> l1 = new ArrayList<Integer>(); // make a copy
        l1.addAll(l);

        ArrayList<Float>   w1 = new ArrayList<Float>(); // make a copy
        w1.addAll(w);

        // offsets for the suppressed locations after random sampling
//        int[][] nbhood = getoffsetsxyz(20, 5);

        Stepper stpr = new Stepper(radius, is2d, kappa, -1);

        Xk.clear();

        int i = 0;
        while (i < n && w1.size()>0) {

            ArrayList<Double> c1 = new ArrayList<Double>(l1.size());

            c1.clear();
            for (int j = 0; j < w1.size(); j++) c1.add((double) w1.get(j) + ((j > 0) ? c1.get(j - 1) : 0));

            // sample one
            double totalmass = c1.get(c1.size()-1);
            double u1 = (totalmass) * new Random().nextDouble(); // / (float)n , no division with n

            int k = 0;
            while (u1 > c1.get(k) && k<c1.size()-1) k++;

            // add it to the output list
            int ii = l1.get(k);
            int x = ii % N;
            int z = ii / (N * M);
            int y = ii / N - z * M;

            float tmin = Float.POSITIVE_INFINITY;
            float tmax = Float.NEGATIVE_INFINITY;

            Arrays.fill(stpr._pcws0, 0);

            for (int j = 0; j < stpr._pcws0.length; j++) { // as used in _iter0 and _iter1

                int xi = x + stpr._p0[j][0];
                int yi = y + stpr._p0[j][1];
                int zi = z + stpr._p0[j][2];
                int li = zi*(N*M)+yi*N+xi;
                boolean inimg = xi>=0 && xi<N && yi>=0 && yi<M && zi>=0 && zi<P;

                if (inimg) { // remove sampled particles from l1,w1, so that the next sampling goes elsewhere
                    int remidx = l1.indexOf(li);
                    if (remidx!=-1) {
                        l1.remove(remidx); // won't be chosen as seed in this round
                        w1.remove(remidx);
                    }
                }

                stpr._pcws0[j] = (inimg && suppmap[li]==0)?tness[li]:0; // suppmap[li]==0 only here so that each object gets roughly ro particles

                if (stpr._pcws0[j]<tmin) tmin = stpr._pcws0[j];
                if (stpr._pcws0[j]>tmax) tmax = stpr._pcws0[j];

            }

            for (int j = 0; j < stpr._pcws0.length; j++) {
                stpr._pcws0[j] = (tmax-tmin>Float.MIN_VALUE)?((stpr._pcws0[j]-tmin)/(tmax-tmin)):0;
                stpr._pcws0[j] = (float) (Math.pow(stpr._pcws0[j],weight_deg) * 1f);
                stpr._pcws0[j] = stpr._pcws0[j] + ((j==0)? 0 : stpr._pcws0[j-1]);
            }

            float wmass = stpr._pcws0[stpr._pcws0.length-1];
            if (wmass<=Float.MIN_VALUE) {
                IJ.log("wmass["+x+","+y+","+z+";"+radius+"]=0");
                suppmap[ii] = -1;
                continue;
            }

            // pick ro values and add to Xk based on _pcws0 importance sampling
            float uu_1 = (wmass/ro)  * rndgen.nextFloat();
            int s = 0;

            int prevsize = Xk.size();
            float sumweights = 0;
            float xc = 0, yc = 0, zc=0;

            for (int j = 0; j < ro; j++) { // ro samples per object

                float uu_j = uu_1 + j * (wmass/ro);

                while (uu_j>stpr._pcws0[s] && s<stpr._pcws0.length-1) s++;

                // add s-th
                int xj = x + stpr._p0[s][0];
                int yj = y + stpr._p0[s][1];
                int zj = z + stpr._p0[s][2];
                int lj = zj*(N*M)+yj*N+xj;
                boolean inimg = xj>=0 && xj<N && yj>=0 && yj<M && zj>=0 && zj<P;

                if (inimg && suppmap[lj]==0) {

                    X particle = new X(xj,yj,zj, Float.NaN,Float.NaN,Float.NaN, tness[lj], tness[lj]);
                    particle.tag = 0;
                    Xk.add(particle);

                    sumweights += tness[lj];
                    xc += tness[lj] * xj;
                    yc += tness[lj] * yj;
                    zc += tness[lj] * zj;

                }
                else {
                    IJ.log("SAMPLE out of img or on the supressed voxel.");
                }

            }

            if (sumweights<=Float.MIN_VALUE) {IJ.log("sumweights<=Float.MIN_VALUE");}

            xc /= sumweights;
            yc /= sumweights;
            zc /= sumweights;

            for (int j = Xk.size()-1; j >= prevsize; j--) {
                // set directions wrt centroid
                Xk.get(j).vx = Xk.get(j).x - xc;
                Xk.get(j).vy = Xk.get(j).y - yc;
                Xk.get(j).vz = Xk.get(j).z - zc;
                float dd = (float) Math.sqrt(Math.pow(Xk.get(j).vx,2)+ Math.pow(Xk.get(j).vy,2)+ Math.pow(Xk.get(j).vz,2));
                Xk.get(j).vx = Xk.get(j).vx / dd;
                Xk.get(j).vy = Xk.get(j).vy / dd;
                Xk.get(j).vz = Xk.get(j).vz / dd;

                if (Float.isNaN(Xk.get(j).vx) || Float.isNaN(Xk.get(j).vy) || Float.isNaN(Xk.get(j).vz)) {
                    IJ.log("at initialization");
                    IJ.log("Float.isNaN(Xk.get(j).vx) || Float.isNaN(Xk.get(j).vy) || Float.isNaN(Xk.get(j).vz)");
                    IJ.log("Xk location= "+Xk.get(j).x+","+Xk.get(j).y+","+Xk.get(j).z);
                    IJ.log("centroid= " +xc+","+yc+","+zc);
                }

                // normalize weights (~tubularity, and sum up to 1 for each object)
                Xk.get(j).w = Xk.get(j).w / sumweights;
            }

            N_o.add(new int[]{x, y, z});

            i++;

        }

        String log = "_init,\t";
        log+="|X|="+Xk.size()+", ";
        int prevYsize = Y.size();
        // at iter0 Xk are used for the estimation (clustering)
        phdmass = estimate(Xk, kernel_radius, suppmap, N, M, P, Y);

        log+="PHD=" + IJ.d2s(phdmass,2)+", ";
        log+="|Yk|=" + (Y.size()-prevYsize) + ", ";
        npcles = Math.round(phdmass)*ro; // will be used to resample in next iteration
        log+="npcles=" + npcles;
        if (verbose) IJ.log(log);

    }

    public boolean _iter1(int N, int M, int P, float[] tness, int[] suppmap) {

        String eventlog = "|X|="+ Xk.size() +" ";

        //-----------------------------------------------------------------------------
        //*** prediction ***//
        int xi, yi, zi, xj, yj, zj, ii;
        float tnessval;
        X particle;

        XPk.clear();
        XPk_cws.clear();
        ZPk.clear();

        ArrayList<X> zpick = new ArrayList<X>();

        for (int i = 0; i < Xk.size(); i++) {

            X xp = Xk.get(i);

            if (Float.isNaN(xp.vx) || Float.isNaN(xp.vy) || Float.isNaN(xp.vz)) {IJ.log("Float.isNaN(xp.vx) || Float.isNaN(xp.vy) || Float.isNaN(xp.vz)=TRUE");}

            xi = Math.round(xp.x);
            yi = Math.round(xp.y);
            zi = Math.round(xp.z);

            int ip = mm.getdirection(xp.vx, xp.vy, xp.vz); // direction

            //-----------------------------------------------------------------------------
            Arrays.fill(mm.pcws, 0); // ZPk

            float tmin = Float.POSITIVE_INFINITY;
            float tmax = Float.NEGATIVE_INFINITY;

            for (int j = 0; j < mm.sz; j++) {
                xj = xi + mm.p[j][0];
                yj = yi + mm.p[j][1];
                zj = zi + mm.p[j][2];
                ii = zj*(N*M)+yj*N+xj;

                mm.pcws[j] = (xj>=0 && xj<N && yj>=0 && yj<M && zj>=0 && zj<P)? tness[ii] : 0; // && suppmap[ii]==0

                if (mm.pcws[j]<tmin) tmin = mm.pcws[j];
                if (mm.pcws[j]>tmax) tmax = mm.pcws[j];
            }

            float tsum = 0;
            for (int j = 0; j < mm.sz; j++) {
                mm.pcws[j] = (tmax-tmin> Float.MIN_VALUE)? ((mm.pcws[j]-tmin)/(tmax-tmin)) : (1f/mm.sz);
                mm.pcws[j] = (float) (Math.pow(mm.pcws[j],weight_deg) * mm.w[ip][j]);
                tsum += mm.pcws[j];
            }

            if (tsum<=Float.MIN_VALUE) continue; // go to the next Xk particle predictions

            for (int j = 0; j < mm.sz; j++) {
                mm.pcws[j] = (mm.pcws[j]/tsum) + ((j==0)? 0 : mm.pcws[j-1]);
            }

            //** ZPk particles **//
            float wmass = mm.pcws[mm.sz-1];
//            IJ.log("check wmass 2 : " + wmass);
            float u1 = (wmass/ni)  * rndgen.nextFloat();

            int s = 0;

            zpick.clear();

            for (int j = 0; j < ni; j++) { // predict ni ZPk from each Xk particle

                float uj = u1 + j * (wmass/ni);

                while (uj>mm.pcws[s] && s<mm.sz-1) s++;

                xj = xi + mm.p[s][0]; // add s-th
                yj = yi + mm.p[s][1];
                zj = zi + mm.p[s][2];

                if (xj>=0 && xj<N && yj>=0 && yj<M && zj>=0 && zj<P) {

                    ii = zj*(N*M)+yj*N+xj;

                    tnessval = tness[ii];
                    // finally select those that were not suppressed (that don't overlap with previous particle trace)
//                    if (suppmap[ii]==0) {
                    particle = new X(xj,yj,zj, mm.u[s][0],mm.u[s][1],mm.u[s][2], 1f, tnessval);
                    particle.tag = xp.tag;
                    zpick.add(particle);
//                    }
                }
            }

            // add once it is clear how many there are
            for (int j = 0; j < zpick.size(); j++) {
                ZPk.add(zpick.get(j));
            }

            //-----------------------------------------------------------------------------
            Arrays.fill(mm.pcws, 0); // XPk

            tmin = Float.POSITIVE_INFINITY;
            tmax = Float.NEGATIVE_INFINITY;

            for (int j = 0; j < mm.sz; j++) {
                xj = xi + mm.p[j][0];
                yj = yi + mm.p[j][1];
                zj = zi + mm.p[j][2];
                ii = zj*(N*M)+yj*N+xj;

                mm.pcws[j] = (xj>=0 && xj<N && yj>=0 && yj<M && zj>=0 && zj<P)? tness[ii] : 0; // && suppmap[ii]==0

                if (mm.pcws[j]<tmin) tmin = mm.pcws[j];
                if (mm.pcws[j]>tmax) tmax = mm.pcws[j];
            }

            tsum = 0;
            for (int j = 0; j < mm.sz; j++) {
                mm.pcws[j] = (tmax-tmin> Float.MIN_VALUE)? ((mm.pcws[j]-tmin)/(tmax-tmin)) : (1f/mm.sz);
                mm.pcws[j] = (float) (Math.pow(mm.pcws[j],weight_deg) * mm.w[ip][j]);
                tsum += mm.pcws[j];
            }

            if (tsum<=Float.MIN_VALUE) continue; // go to the next Xk particle

            for (int j = 0; j < mm.sz; j++) {
                mm.pcws[j] = (mm.pcws[j]/tsum) + ((j==0)? 0 : mm.pcws[j-1]);
            }

            //** XPk particles **//
            wmass = mm.pcws[mm.sz-1];
//            IJ.log("check wmass 1 : " + wmass);
            u1 = (wmass/ni)  * rndgen.nextFloat();

            s = 0;

            ArrayList<X> xpick = new ArrayList<X>();
            xpick.clear();

            for (int j = 0; j < ni; j++) { // predict ni XPk from each Xk particle

                float uj = u1 + j * (wmass/ni);

                while (uj>mm.pcws[s] && s<mm.sz-1) s++;

                xj = xi + mm.p[s][0]; // add s-th
                yj = yi + mm.p[s][1];
                zj = zi + mm.p[s][2];

                if (xj>=0 && xj<N && yj>=0 && yj<M && zj>=0 && zj<P) {

                    ii = zj*(N*M)+yj*N+xj;

                    tnessval = tness[ii];
                    // finally select those that were not suppressed (that don't overlap with previous particle trace)
//                    if (suppmap[ii]==0) {
                    particle = new X(xj,yj,zj, mm.u[s][0],mm.u[s][1],mm.u[s][2], 1f, tnessval);
                    particle.tag = xp.tag;
                    xpick.add(particle);
//                    }
                }
            }

            // add once it is clear how many there are
            for (int j = 0; j < xpick.size(); j++) {
                xpick.get(j).w = pS * xp.w * (1f/xpick.size());
                if (Float.isNaN(xpick.get(j).vx) || Float.isNaN(xpick.get(j).vy) || Float.isNaN(xpick.get(j).vz))
                {IJ.log("Float.isNaN(xpick.get(j).vx) || Float.isNaN(xpick.get(j).vy) || Float.isNaN(xpick.get(j).vz)=TRUE");}
                XPk.add(xpick.get(j));
                if (XPk_cws.size()==0)  XPk_cws.add(xpick.get(j).w);
                else                    XPk_cws.add(xpick.get(j).w + XPk_cws.get(XPk_cws.size()-1));
            }

        } // go through Xk PHD particles

        eventlog += "|XPk|= "+XPk.size()+" ";//+ "  XPk_cws[end]="+XPk_cws.get(XPk_cws.size()-1)+"  ZPk.size()="+ZPk.size()+"  ";

        if (XPk.size()==0) {IJ.log("XPk.size()==0"); return false;}
        if (ZPk.size()==0) {IJ.log("ZPk.size()==0"); return false;}

        //-----------------------------------------------------------------------------
        //** measure **//
        Zk.clear();
        measure(ZPk, kernel_radius, Y, N, M, tness, suppmap, Zk);
        eventlog += "|Z|=" + Zk.size() + " ";

        if (Zk.size()==0) {IJ.log("Zk.size()==0"); return false;}

        //** update **//
        XPk_cws = update(XPk, Zk); // updated weights will be added to cws output

        phdmass = XPk_cws.get(XPk_cws.size()-1); eventlog += "PHD="  +IJ.d2s(phdmass,2) + " ";

        if (Math.round(phdmass)<1) {IJ.log("Math.round(phdmass)<1"); return false;}
        if (Math.round(phdmass)>OBJECT_LIMIT) {IJ.log("#OBJ. LIMIT : Math.round(phdmass)>" + OBJECT_LIMIT); phdmass = OBJECT_LIMIT;}

        int prevYsize = Y.size();

        //** estimate **//
        ArrayList<X> Rk = new ArrayList<X>();
        ArrayList<Float> Rcws = new ArrayList<Float>();  // XPk subsets used in the estimation
        est(XPk, kernel_radius, MIN_CLUST_SIZE, Integer.MAX_VALUE, wmin, N, M, suppmap, Y, Rk, Rcws);

        if (Rk.size()==0) {IJ.log("Rk.size()==0"); return false;}

        eventlog += "|Y|=" + (Y.size()-prevYsize)+" ";

        //** resample **//
        npcles = Math.round(phdmass)*ro;
        resample(Rk, Rcws, npcles, phdmass/(float)npcles, Xk);

        // now fill out the suppmap[] knowing XPk, Rk and Xk in the same way it was filled in est(),
        // est() needs one more output ArrayList<Integer> with Y node indexes of each Rk
        // apply current XPk to the suppmap

        if (verbose) IJ.log(eventlog);

        return true;

    }

    private void measure(
            ArrayList<X> zp,
            float dist,
            ArrayList<Node> Yprev,
            int N, int M, //int P,
            float[] tness,
            int[] suppmap,
            ArrayList<X> zout//,
            //ArrayList<ArrayList<Integer>> zc
    ) {

        float gauss_limit = 3*dist;
        Arrays.fill(xw, 0);

        for (int i = 0; i < zp.size(); i++) {

            for (int j = i; j < zp.size(); j++) {

                if (i==j) {

                    xw[i] += zp.get(i).tness; // use tness instead of weight x.get(i).w

                }
                else {

                    float d2 = (float) Math.pow(zp.get(i).x - zp.get(j).x, 2);

                    if (d2 < Math.pow(gauss_limit,2)) {

                        d2 += Math.pow(zp.get(i).y - zp.get(j).y,2);

                        if (d2 < Math.pow(gauss_limit,2)) {

                            d2 += Math.pow(zp.get(i).z - zp.get(j).z,2);

                            if (d2 < Math.pow(gauss_limit,2)) {

                                float dcost = (float) Math.exp(-d2/(2*Math.pow(dist,2)));
                                xw[i] += zp.get(j).tness * dcost; // x.get(j).w
                                xw[j] += zp.get(i).tness * dcost; // x.get(i).w

                            }

                        }

                    }

                }

            }

        }

        meanShift(zp, xw, dist);     // conv filled up
        clustering(zp.size(), 2f);   // labels filled up, nbridxs filled up
        group_measurement(zp, zp.size(), MIN_CLUST_SIZE, Integer.MAX_VALUE, Yprev, N, M, tness, suppmap, zout); // zc

    }

    private void group_measurement(
            ArrayList<X> zp,
            int zpLen,
            int nclust_min,
            int nclust_max,
            ArrayList<Node> Yprev,
            int N, int M,
            float[] tness,
            int[] suppmap,
            ArrayList<X> z//,
//            ArrayList<ArrayList<Integer>> zc
    ) {

        // no need for estimation of the centroid, conv[][] is not used, meanShift() has clustering role
        Arrays.fill(checked, false);

        for (int i = 0; i < ZPCk.size(); i++) ZPCk.get(i).clear();
        ZPCk.clear();

        ArrayList<Integer> cnt = new ArrayList<Integer>(); // only used for sorting

        for (int i = 0; i < zpLen; i++) {

            if (!checked[i]) {
                checked[i] = true;

                ArrayList<Integer> clsidxs = new ArrayList<Integer>();
                clsidxs.add(i);
                boolean notSupp = suppmap[Math.round(zp.get(i).z)*(N*M)+Math.round(zp.get(i).y)*N+Math.round(zp.get(i).x)]==0;

                for (int j = i+1; j < zpLen; j++) { // check the rest
                    if (!checked[j] && labels[j]==labels[i]) {
                        checked[j] = true;
                        clsidxs.add(j);
                        notSupp = notSupp || (suppmap[Math.round(zp.get(j).z)*(N*M)+Math.round(zp.get(j).y)*N+Math.round(zp.get(j).x)]==0);
                    }
                }

                if (clsidxs.size()>=nclust_min  && notSupp) { // cluster had enough samples and at least one that did not overlap with suppmap[]
                    cnt.add(clsidxs.size());
                    ZPCk.add(clsidxs);
                }
            }
        }

        // sort clusters by the amount of particles in the cluster
        int[] desc_idx = Tools.descending(cnt); // cnt will be modified too todo: not necessary really

        ArrayList<Float> zsc        = new ArrayList<Float>(ZPCk.size());        // score per cluster
        ArrayList<Integer> zsc_idx  = new ArrayList<Integer>(ZPCk.size());  // selected index per cluster

        for (int ii=0; ii<Math.min(nclust_max,ZPCk.size()); ii++) {

            int ci = desc_idx[ii]; // go from the highest, take cluster index

            float sc = Float.POSITIVE_INFINITY;
            int sc_idx = -1;

            for (int i = 0; i < ZPCk.get(ci).size(); i++) { // go through elements of one cluster, there should be at least one that's not overlapping with suppmap[]

                int xii = ZPCk.get(ci).get(i);
                X xi = zp.get(xii);

                if (suppmap[Math.round(xi.z)*(N*M)+Math.round(xi.y)*N+Math.round(xi.x)]==0) { // should happen at least once

                    float a1x = Yprev.get(xi.tag).loc[0];
                    float a1y = Yprev.get(xi.tag).loc[1];
                    float a1z = Yprev.get(xi.tag).loc[2];

                    float a2x = xi.x;
                    float a2y = xi.y;
                    float a2z = xi.z;

                    float score = 0;

                    for (int j = 0; j < ZPCk.get(ci).size(); j++) {
                        if (j!=i) {

                            X xj = zp.get(ZPCk.get(ci).get(j));

                            float a0x = xj.x;
                            float a0y = xj.y;
                            float a0z = xj.z;

                            score += d2(a1x, a1y, a1z, a2x, a2y, a2z, a0x, a0y, a0z);

                        }
                    }

                    if (score<sc) {
                        sc = score;   // update score
                        sc_idx = xii; // update index
                    }
                } // if it was not overlapping
            } // loop cluster

            zsc.add(sc);
            zsc_idx.add(sc_idx);

        } // loop clusters

        if (zsc.size()==0) { // stop without extracting z
            return;
        }

        // store the sampled observations into z, those that had clusters with at least one not on the suppmap (XPk subset used to get Y fills the suppmap)
        z.clear();

        for (int i = 0; i < zsc_idx.size(); i++) {

            int x_Z = Math.round(zp.get( zsc_idx.get(i) ).x);
            int y_Z = Math.round(zp.get( zsc_idx.get(i) ).y);
            int z_Z = Math.round(zp.get( zsc_idx.get(i) ).z);

            X xtt = new X(zp.get(zsc_idx.get(i)));
            xtt.tness = tness[z_Z*(N*M)+y_Z*(N)+x_Z];
            z.add(xtt);

//            if (x_Z>=0 && x_Z<N && y_Z>=0 && y_Z<M && z_Z>=0 && z_Z<P)
//                if (suppmap[z_Z*(N*M)+y_Z*(N)+x_Z]==0 || true)
        }

    }

    private float d2(
            float a1x, float a1y, float a1z,
            float a2x, float a2y, float a2z,
            float a0x, float a0y, float a0z
            ) {

        float a21x = a2x-a1x;
        float a21y = a2y-a1y;
        float a21z = a2z-a1z;

        float a01x = a0x-a1x;
        float a01y = a0y-a1y;
        float a01z = a0z-a1z;

        float a02x = a0x-a2x;
        float a02y = a0y-a2y;
        float a02z = a0z-a2z;

        // d2=0 if it is outside
        if (a21x*a01x+a21y*a01y+a21z*a01z<=0)
            return (a01x*a01x+a01y*a01y+a01z*a01z);//Float.POSITIVE_INFINITY;
        if ((-a21x)*a02x+(-a21y)*a02y+(-a21z)*a02z<=0)
            return (a02x*a02x+a02y*a02y+a02z*a02z);//Float.POSITIVE_INFINITY;

//        return 1;

        // point to line distance 3D
        // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

        float a10x = a1x-a0x;
        float a10y = a1y-a0y;
        float a10z = a1z-a0z;

        float[] a21_x_a10 = new float[3]; // this can be a global variable

        crossprod(a21x, a21y, a21z, a10x, a10y, a10z, a21_x_a10);

        return (float) (
                (Math.pow(a21_x_a10[0],2)+Math.pow(a21_x_a10[1],2)+Math.pow(a21_x_a10[2],2)) / ((a21x*a21x)+(a21y*a21y)+(a21z*a21z))
        );

    }

    private void crossprod(float a1, float a2, float a3, float b1, float b2, float b3, float[] v) {
        // v is cross product of (a1, a2, a3) and (b1, b2, b3)
        v[0] = a2*b3 - b2*a3;
        v[1] = -(a1*b3-b1*a3);
        v[2] = a1*b2-b1*a2;
    }

    private void estimate(
            ArrayList<ArrayList<Integer>> xg,       // clusters of the particle indexes
            ArrayList<X> x,                         //
            int N,
            int M,
            int P,
            int[] suppmap,
            ArrayList<X> Gk,                        // output select particles that are clustered (in xg) and not suppressed
            ArrayList<Float> Gcws,                  // output cumulative weight sum
            ArrayList<Node> nout                    // output node list with estimations
    ) {

        Gk.clear();
        Gcws.clear();

        for (int i = 0; i < xg.size(); i++) {

            float wsum = 0, cx=0, cy=0, cz=0;
            ArrayList<Integer> tags = new ArrayList<Integer>();

            for (int j = 0; j < xg.get(i).size(); j++) {

                int xii = xg.get(i).get(j);

                wsum += x.get(xii).w;
                cx += x.get(xii).x * x.get(xii).w;
                cy += x.get(xii).y * x.get(xii).w;
                cz += x.get(xii).z * x.get(xii).w;

                tags.add(x.get(xii).tag);

            }

            if (wsum>Float.MIN_VALUE) { // updated weights of the cluster should not sum up to 0

                cx /= wsum;
                cy /= wsum;
                cz /= wsum;

                // remove duplicate tags
                Set<Integer> set = new HashSet<Integer>();
                set.addAll(tags);
                tags.clear();
                tags.addAll(set);

                Node nn = new Node(cx, cy, cz, 1f); // rr
                int newtag = nout.size();
                nout.add(nn);

                for (int j = 0; j < tags.size(); j++) {
                    nout.get(newtag).nbr.add(tags.get(j));
                    nout.get(tags.get(j)).nbr.add(newtag);
                }

                // fill out the suppression map, go through the i-th cluster particle-indexes again
                for (int j = 0; j < xg.get(i).size(); j++) {

                    // assign new tag to the particles of the clusters from XPk
                    x.get(xg.get(i).get(j)).tag = newtag;// XPk particle will get a tag that corresponds to the estimate

                    int idx = xg.get(i).get(j);
                    int x1 = Math.round(x.get(idx).x);
                    int y1 = Math.round(x.get(idx).y);
                    int z1 = Math.round(x.get(idx).z);

                    if (x1>=0 && x1<N && y1>=0 && y1<M && z1>=0 && z1<P) {
                        int idx1 = z1*(N*M)+y1*N+x1;
                        if (suppmap[idx1]==0) {

                            // *** add Gk and Gcws
                            Gk.add(x.get(idx));

                            if (Gcws.size()==0) Gcws.add(x.get(idx).w);
                            else Gcws.add(x.get(idx).w + Gcws.get(Gcws.size() - 1));

                            // ***
                            suppmap[idx1] = newtag;

                        }
                    }





                }



            }
            else {
                IJ.log("cluster " + i + " was removed from xg, due to weight sum zero.");
                xg.get(i).clear(); // cancel those that had weight sum 0 so that they are not used later, they won't get a new tag either
            }

        }

        // cluster weighted phd particles using mean-shift and get the weighted mean out of each cluster
//        float outcws = 0;
//        for (int i = 0; i < x.size(); i++) {
//            xw[i] = x.get(i).tness; // used to be w
//            outcws += x.get(i).w;
//        }
//        meanShift(x, xw, dist);
//        clustering(x.size(), 2f);
//        group(x, MIN_CLUST, Integer.MAX_VALUE, suppmap, rsupp, N, M, P, y); // estimate() uses x phd particles
//        return outcws;

    }

    private void est(
            ArrayList<X> x,
            float kradius,
            int nclust_min,
            int nclust_max,
            float wmin,
            int N,
            int M,
            int[] suppmap,
            ArrayList<Node> y,
            ArrayList<X> Rk,
            ArrayList<Float> Rcws
    ) {

        // cluster weighted phd particles
        for (int i = 0; i < x.size(); i++) {xw[i] = x.get(i).w;}

        meanShift(x, xw, kradius);
        clustering(x.size(), 2f);
        group_estimate(x, nclust_min, nclust_max, wmin, N, M, suppmap, y, Rk, Rcws);

    }

    private float estimate(ArrayList<X> x, float dist, int[] suppmap, int N, int M, int P, ArrayList<Node> y) {

        // cluster weighted phd particles using mean-shift and get the weighted mean out of each cluster
        float outcws = 0;
        for (int i = 0; i < x.size(); i++) {
            xw[i] = x.get(i).tness; // used to be w
            outcws += x.get(i).w;
        }

        meanShift(x, xw, dist);
        clustering(x.size(), 2f);
        group(x, 0, Integer.MAX_VALUE, suppmap, N, M, P, y); // estimate() uses x phd particles   rsupp

        return outcws;
    }

    private void clustering(int xlen, float dist) {

        // indxs represent indexes of values that need to be clustered according to their values read before
        // threshold_dists is the distance limit
        // output is list of unique labels
//        int[] labels = new int[values.length];

        float dist2 = dist*dist;

        for (int i = 0; i < xlen; i++) {
            labels[i] = i;
            nbridxs[i].clear();
        }

        for (int i = 0; i < xlen; i++) {

            for (int j = i+1; j < xlen; j++) {

                float d2 = (float) Math.pow(conv[i][0]-conv[j][0],2);

                if (d2<dist2) {

                    d2 += Math.pow(conv[i][1]-conv[j][1],2);

                    if (d2<dist2) {

                        d2 += Math.pow(conv[i][2]-conv[j][2],2);

                        if (d2<dist2) {

                            nbridxs[i].add(j);
                            nbridxs[j].add(i);

                        }

                    }

                }

            }

//            for (int j = i; j < xlen; j++) {
//                if(i==j) {
//                    dists[i][j] = 0;
//                }
//                else {
//                    float x2 = conv[i][0]-conv[j][0]; x2 = x2*x2;
//                    float y2 = conv[i][1]-conv[j][1]; y2 = y2*y2;
//                    float z2 = conv[i][2]-conv[j][2]; z2 = z2*z2;
//                    dists[i][j] = x2+y2+z2;
//                    dists[j][i] = dists[i][j];
//                }
//            }

        }

        for (int i = 0; i < xlen; i++) {

            for (int nbri = 0; nbri < nbridxs[i].size(); nbri++) {

                int j = (Integer) nbridxs[i].get(nbri);

                // propagate labels
                if (labels[j] != labels[i]) {

                    int currLabel = labels[j];
                    int newLabel  = labels[i];

                    labels[j] = newLabel;

                    //set all that also were currLabel to newLabel
                    for (int k = 0; k < xlen; k++)
                        if (labels[k]==currLabel)
                            labels[k] = newLabel;

                }

            }



//            // one versus the rest
//            for (int j = 0; j < xlen; j++) {
//                if (i != j) {
//                    if (dists[i][j]<=dist2) {
//
//
//
//
//                        // exchange labels
//                        if (labels[j] != labels[i]) {
//
//                            int currLabel = labels[j];
//                            int newLabel  = labels[i];
//
//                            labels[j] = newLabel;
//
//                            //set all that also were currLabel to newLabel
//                            for (int k = 0; k < xlen; k++)
//                                if (labels[k]==currLabel)
//                                    labels[k] = newLabel;
//
//                        }
//
//                    }
//                }
//            }




        }

//        return labels;

    }

    private void meanShift(ArrayList<X> xin, float[] xw, float dist) {

//        float[][] conv = new float[xin.size()][3];

        for (int i = 0; i < xin.size(); i++) {
            conv[i][0] = xin.get(i).x;
            conv[i][1] = xin.get(i).y;
            conv[i][2] = xin.get(i).z;
        }

        float[] new_v = new float[3];

        for (int i = 0; i < xin.size(); i++) {

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

//        return conv;

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

    private void group(ArrayList<X> xin, int count_min, int nclust_max, int[] suppmap, int N, int M, int P, ArrayList<Node> nout) { // int rsupp,

        Arrays.fill(checked, false);

        float wsum, cx, cy, cz;
        ArrayList<Integer> tags = new ArrayList<Integer>();
        ArrayList<Integer> xidx = new ArrayList<Integer>();
        int nr_clusters = 0;

        for (int i = 0; i < xin.size(); i++) {

            if (!checked[i]) {

                wsum = xin.get(i).tness; // w
                cx = xin.get(i).x * xin.get(i).tness; // w
                cy = xin.get(i).y * xin.get(i).tness; // w
                cz = xin.get(i).z * xin.get(i).tness; // w

                nr_clusters++;

                xidx.clear();
                xidx.add(i);

                tags.clear();
                if (xin.get(i).tag>0)
                    tags.add(xin.get(i).tag);

                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < xin.size(); j++) {
                    if (!checked[j] && labels[j]==labels[i]) {

                        wsum += xin.get(j).tness; // w

                        cx += xin.get(j).x * xin.get(j).tness; // w
                        cy += xin.get(j).y * xin.get(j).tness; // w
                        cz += xin.get(j).z * xin.get(j).tness; // w

                        xidx.add(j);

                        if (xin.get(j).tag>0)
                            tags.add(xin.get(j).tag);

                        count++;
                        checked[j] = true;

                    }
                }

//                if (wsum>Float.MIN_VALUE) {
                    if (count >= count_min) {

                        if (wsum>Float.MIN_VALUE) {

                        cx /= wsum;
                        cy /= wsum;
                        cz /= wsum;
//                        float rr = 0;
//                        for (int j = 0; j < xidx.size(); j++) {
//                            float px = xin.get(xidx.get(j)).x;
//                            float py = xin.get(xidx.get(j)).y;
//                            float pz = xin.get(xidx.get(j)).z;
//                            float pw = xin.get(xidx.get(j)).w;
//                            rr += Math.sqrt(Math.pow(px-cx,2) + Math.pow(py-cy,2) + Math.pow(pz-cz,2)) * pw;
//                        }

//                        rr /= wsum;
                        //rr = Math.round(rr); // not necessary as it;s not used to fill any pixel map
                        //rr = (rr<1)?1:rr;
                        //rr = (rr> MultiTT.MAX_RADIUS)?MultiTT.MAX_RADIUS:rr;

                        // remove duplicate tags
                        Set<Integer> set = new HashSet<Integer>();
                        set.addAll(tags);
                        tags.clear();
                        tags.addAll(set);

                        Node nn = new Node(cx, cy, cz, 2f); // rr
                        int newtag = nout.size();
                        nout.add(nn);

                        for (int j = 0; j < tags.size(); j++) {
                            nout.get(newtag).nbr.add(tags.get(j));
                            nout.get(tags.get(j)).nbr.add(newtag);
                        }

                            for (int j = 0; j < xidx.size(); j++) {
                            // assign new tag to the particles of the cluster
                            xin.get(xidx.get(j)).tag = newtag;

                                for (int k = 0; k < offxyz[R_supp].length; k++) {

                                    int xi = Math.round(xin.get(xidx.get(j)).x) + offxyz[R_supp][k][0];
                                    int yi = Math.round(xin.get(xidx.get(j)).y) + offxyz[R_supp][k][1];
                                    int zi = Math.round(xin.get(xidx.get(j)).z) + offxyz[R_supp][k][2];

                                    if (xi>=0 && xi<N && yi>=0 && yi<M && zi>=0 && zi<P) {
                                        int ii = zi*(N*M)+yi*N+xi;
                                        if (suppmap[ii]==0)
                                            suppmap[ii] = newtag;
                                    }
                                }
                            }


                        }

                    } // count >= count_min
//                    else {
//                        IJ.log("cluster with  "+xidx.size()+" elements set w and ta to zero");
//                        // not enough counts, set weights to 0 so that they're surely not resampled later on
//                        for (int j = 0; j < xidx.size(); j++) {
//                            xin.get(xidx.get(j)).w = 0;
//                            xin.get(xidx.get(j)).tag = 0;
//                        }
//
//                    }

//                }
//                else {
                    //IJ.log("there was cluster with enough counts but zero wsum!"); // won't be resampled
//                }
            }
        }
    }

    private void group_estimate(
            ArrayList<X> xin,
            int count_min,
            int nclust_max,
            float wmin,
            int N,
            int M,
            int[] suppmap,
            ArrayList<Node> nout,
            ArrayList<X> Rk,
            ArrayList<Float> Rcws
    ) {

        Rk.clear();
        Rcws.clear();

        Arrays.fill(checked, false);

        float wsum, cx, cy, cz;
        ArrayList<Integer> tags = new ArrayList<Integer>();
        ArrayList<Integer> xidx = new ArrayList<Integer>();
        int nr_clusters = 0;

        for (int i = 0; i < xin.size(); i++) {

//            suppmap[Math.round(xin.get(i).z)*(N*M)+Math.round(xin.get(i).y)*N+Math.round(xin.get(i).x)] = xin.get(i).tag;

            if (!checked[i]) {

                wsum = xin.get(i).w;
                cx = xin.get(i).x * xin.get(i).w; // tness?
                cy = xin.get(i).y * xin.get(i).w;
                cz = xin.get(i).z * xin.get(i).w;

                nr_clusters++;

                xidx.clear();
                xidx.add(i);

                tags.clear();
                if (xin.get(i).tag>0)
                    tags.add(xin.get(i).tag);

                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < xin.size(); j++) {
                    if (!checked[j] && labels[j]==labels[i]) {

                        wsum += xin.get(j).w;

                        cx += xin.get(j).x * xin.get(j).w;
                        cy += xin.get(j).y * xin.get(j).w;
                        cz += xin.get(j).z * xin.get(j).w;

                        xidx.add(j);

                        if (xin.get(j).tag>0)
                            tags.add(xin.get(j).tag);

                        count++;
                        checked[j] = true;

                    }
                }

                if (count>=count_min) {

                    if (wsum>=wmin) {

                        cx /= wsum;
                        cy /= wsum;
                        cz /= wsum;

                        Set<Integer> set = new HashSet<Integer>();
                        set.addAll(tags);
                        tags.clear();
                        tags.addAll(set);

                        //** add the node **//
                        Node nn = new Node(cx, cy, cz, 1f); // rr
                        int newtag = nout.size();
                        nout.add(nn);

                        for (int j = 0; j < tags.size(); j++) {
                            nout.get(newtag).nbr.add(tags.get(j));
                            nout.get(tags.get(j)).nbr.add(newtag);
                        }

                        for (int j = 0; j < xidx.size(); j++) {

                            //** assign new tag to the particles of the cluster **//
                            xin.get(xidx.get(j)).tag = newtag;

                            if (Float.isNaN(xin.get(xidx.get(j)).vx) || Float.isNaN(xin.get(xidx.get(j)).vy) || Float.isNaN(xin.get(xidx.get(j)).vz))
                            {IJ.log("Float.isNaN(xin.get(xidx.get(j)).vx) || Float.isNaN(xin.get(xidx.get(j)).vy) || Float.isNaN(xin.get(xidx.get(j)).vz)=TRUE");}

                            //** add particles used for the estimate (for Xk resampling later) **//
                            Rk.add(xin.get(xidx.get(j)));

                            if (Rcws.size()==0) Rcws.add(xin.get(xidx.get(j)).w);
                            else                Rcws.add(xin.get(xidx.get(j)).w + Rcws.get(Rcws.size() - 1));

                            //** fill suppmap with selected particles **//
                            int xi = Math.round(xin.get(xidx.get(j)).x);
                            int yi = Math.round(xin.get(xidx.get(j)).y);
                            int zi = Math.round(xin.get(xidx.get(j)).z);
                            suppmap[zi*(N*M)+yi*N+xi] = newtag;

                            // for (int k = 0; k < offxyz[R_supp].length; k++)

                        }

                    }

                } // count >= count_min

            }

        }

    }

    private ArrayList<Float> update(ArrayList<X> x, ArrayList<X> z) {// , ArrayList<Integer> cnt // if updated with X instances, outputs CWS with updated weights, for resampling

        g = new float[x.size()][z.size()];
        Cz = new float[z.size()];

        // calculate Cz
        for (int j = 0; j < z.size(); j++) {

            // Cz[j] initialize with clutter PHD here
            Cz[j] = (float) Math.exp(-kclutt*z.get(j).tness); // clutter(z.get(j).tness, kclutt, cluttertness); // *(float)tnessinit

            for (int i = 0; i < x.size(); i++) {

//                if (x.get(i)!=null) { // null are the ones that are discarded at clustering

                    g[i][j] = gzx(z.get(j), x.get(i), gzx_sigma);

                    Cz[j] += pD * g[i][j] * x.get(i).w;

//                }

            }

        }
        
        // update weights using g, Cz and clutter phd
        ArrayList<Float> cws = new ArrayList<Float>(x.size());

        for (int i = 0; i < x.size(); i++) {

//            if (x.get(i)!=null) {

                float wup = (1-pD) * x.get(i).w; // wup == updated weight for i-th particle

                for (int j = 0; j < z.size(); j++) {
                    wup += (pD*g[i][j]*x.get(i).w)/Cz[j]; // update particle weight equation
                }

                x.get(i).w = wup;

                if (cws.size()==0)   cws.add(i, wup);
                else        cws.add(i, wup + cws.get(cws.size()-1));

//            }

        }

        return cws;

    }

    public static int[] cluster(ArrayList<X> Xlist, double ddist) {

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
                                    Math.pow(ddist, 2); // Xlist.get(i).sig+Xlist.get(j).sig


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
//                            cc.sig += vals.get(j).sig;

                            count++;
                            checked[j] = true;

                        }
                    }
                }

                out.add(new Z(cc.x/count, cc.y/count, cc.z/count, count)); // cc.sig/count,
//                out.add(new float[]{centroid/count, count});

            }
        }

        return out;

    }

    private ArrayList<Integer> importsamp(ArrayList<Float> lcws, int n) {
        // systematic resampling, Beyond Kalman Filtering, Ristic et al.
        float totalmass = lcws.get(lcws.size()-1);
        float u1 = (totalmass/(float)n) * rndgen.nextFloat();

        ArrayList<Integer> out = new ArrayList<Integer>(n);
        out.clear();
        int i = 0;
        for (int j = 0; j < n; j++) {
            float uj = u1 + j*(totalmass/(float)n);
            while (uj > lcws.get(i)) i++;
            out.add(i);
        }

        return out;

    }

    private ArrayList<Integer> sequence(int n) {
        ArrayList<Integer> out = new ArrayList<Integer>(n);
        for (int i = 0; i < n; i++) {
            out.add(i);
        }
        return out;
    }

    void resample(ArrayList<X> Xin, ArrayList<Float> Cws, int N, float weight, ArrayList<X> Xout) {
        // systematic resampling algorithm, Beyond Kalman Filtering, Ristic et al.
        float wmass = Cws.get(Cws.size()-1);
        float u1 = (wmass/N) * rndgen.nextFloat();

        Xout.clear();
        
        int i = 0;

        for (int j = 0; j < N; j++) {
            float uj = u1 + j*(wmass/(float)N);
            while (uj > Cws.get(i) && i<Cws.size()-1) i++;
            X tt    = new X(Xin.get(i));
            tt.w    = weight;
            tt.tag  = Xin.get(i).tag;

            if (Float.isNaN(tt.vx) || Float.isNaN(tt.vy) || Float.isNaN(tt.vz))
            {IJ.log("Float.isNaN(tt.vx) || Float.isNaN(tt.vy) || Float.isNaN(tt.vz)=TRUE");}

            Xout.add(tt);
        }

    }

    public static float clutter(float t, float K, float tc) {
        return (float) Math.exp(-K*(t-tc));
    }

    private int get_undiscovered(ArrayList[] discovered){

        for (int i = 0; i < discovered.length; i++) { // first element in stalker list of nodes is nothing
            if (discovered[i]!=null) {
                for (int j = 0; j < discovered[i].size(); j++) {
                    if (!(Boolean) discovered[i].get(j)) {
                        return i;
                    }
                }
            }
        }

        return -1; // all are discovered

    }

    private ArrayList[] init_discovered_list(ArrayList<Node> nlist) {
        // will be used to book-keep discovered traces (graph vertices), all are not discovered at the beginning
        ArrayList[] discovered = new ArrayList[nlist.size()];
        discovered[0] = null;
        for (int i = 1; i < nlist.size(); i++) { // skip the first index as the first element of node list is dummy (indexing starts from 1)
            discovered[i] = new ArrayList<Boolean>(nlist.get(i).nbr.size());
            for (int j = 0; j < nlist.get(i).nbr.size(); j++) {
                discovered[i].add(false);
            }
        }
        return discovered;
    }

    public ArrayList<Node> bfs1(ArrayList<Node> nlist, boolean remove_isolated_tree_with_one_node) {

        /*
        https://en.wikipedia.org/wiki/Breadth-first_search
        1 Breadth-First-Search(Graph, root):
        2
        3     for each node n in Graph:
        4         n.distance = INFINITY
        5         n.parent = NIL
        6
        7     create empty queue Q
        8
        9     root.distance = 0
        10     Q.enqueue(root)
        11
        12     while Q is not empty:
        13
        14         current = Q.dequeue()
        15
        16         for each node n that is adjacent to current:
        17             if n.distance == INFINITY:
        18                 n.distance = current.distance + 1
        19                 n.parent = current
        20                 Q.enqueue(n)
        */

        BfsQueue q = new BfsQueue();

        ArrayList<Node> tree = new ArrayList<Node>();

        int[] dist = new int[nlist.size()];
        Arrays.fill(dist, Integer.MAX_VALUE);
        dist[0] = -1;

        // save indexing in output tree
        int[] nmap = new int[nlist.size()];
        Arrays.fill(nmap, -1);

        // save parent index in current tree
        int[] parent = new int[nlist.size()];
        Arrays.fill(parent,-1);

        tree.add(null);
        int treecnt = 0;

        int seed;

        while ((seed = get_undiscovered(dist))>0) {

            treecnt++;

            dist[seed] = 0;
            nmap[seed] = -1;
            parent[seed] = -1;
            q.enqueue(seed);

            int nodesInTree = 0;

            while (q.hasItems()) {

                // dequeue(), take from FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
                int curr = (Integer) q.dequeue();

                float x = nlist.get(curr).loc[0];
                float y = nlist.get(curr).loc[1];
                float z = nlist.get(curr).loc[2];
                float r = nlist.get(curr).r;
                Node n = new Node(x, y, z, r, (treecnt+1)); // start from RED tree
                if (parent[curr] > 0) n.nbr.add(nmap[parent[curr]]);

                nmap[curr] = tree.size();
                tree.add(n);
                nodesInTree++;

                // for each node adjacent to current
//                int countadj = 0;
                for (int j = 0; j < nlist.get(curr).nbr.size(); j++) {

                    int adj = nlist.get(curr).nbr.get(j);

                    if (dist[adj] == Integer.MAX_VALUE) {
                        dist[adj] = dist[curr] + 1;
                        parent[adj] = curr;
                        // enqueue(), add to FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
                        q.enqueue(adj);
                    }

                }

                // check if there were any neighbours
                if (nodesInTree==1 && !q.hasItems() && remove_isolated_tree_with_one_node) {
                    tree.remove(tree.size()-1);     // remove the one that was just added
                    nmap[curr] = -1;                // cancel the last entry
                }

            }

//            IJ.log("tree[ "+treecnt+" ] ="+nodesInTree+" nodes");
//            if (nodesInTree==1) {
//                tree.remove(tree.size()-1);
//            }

        }

        IJ.log(treecnt+" trees.");

        return tree;

    }

    private int get_undiscovered(int[] dist){

        for (int i = 0; i < dist.length; i++) {
            if (dist[i]>0) {
                if (dist[i]==Integer.MAX_VALUE) {
                    return i;
                }
            }
        }

        return -1;

    }

    public ArrayList<Node> bfs(ArrayList<Node> nlist) {

        /**
         *  breadth-first search (BFS) to traverse the tree from extracted node list
         *  http://en.wikipedia.org/wiki/Breadth-first_search
         *
         1  procedure BFS(G,v) is
         2      let Q be a queue
         3      Q.enqueue(v)
         4      label v as discovered
         5      while Q is not empty
         6         v ← Q.dequeue()
         7         for all edges from v to w in G.adjacentEdges(v) do
         8             if w is not labeled as discovered
         9                 Q.enqueue(w)
         10                label w as discovered
         *
         */

        // will essentially convert ArrayList<Node> with all its linkings (bi-directional connections)
        // into the list of trees ArrayList<ArrayList<Trace>> t where each tree t[i] contains the BFS traverse (with 1 directional connections)
        // BFS needs Queue data structure that will be implemented in class BfsQueue
        // Queue keeps the links between the nodes, link is described as int[] where int[] ~ [curr_node_idx, adjacent_node_idx]
        // knowing just the node index in the queue is not enough, need to know the index of the mother node as well
        // discovered is array of lists that will keep the labels of the discovered adjacent node pairs, bookekeeping for the BFS

        // the key reason for keeping two values is that each time we need the mother index and so for each trace element, including the first node of the trace

        BfsQueue bfsQueue = new BfsQueue();

        ArrayList<Node> tree = new ArrayList<Node>();

        ArrayList[] discovered = init_discovered_list(nlist);

        int[] nodemap = new int[nlist.size()];
        Arrays.fill(nodemap, -1);

        tree.add(null); // size=1
        int tree_count = 1;

        int seed;
        while ((seed = get_undiscovered(discovered))!=-1) {

            float xseed = nlist.get(seed).loc[0];
            float yseed = nlist.get(seed).loc[1];
            float zseed = nlist.get(seed).loc[2];
            float rseed = nlist.get(seed).r;

            nodemap[seed] = tree.size();
            tree.add(new Node(xseed, yseed, zseed, rseed, tree_count));

            // add the neighbors to the queue and label them as discovered
            for (int j = 0; j <nlist.get(seed).nbr.size(); j++) {
                int next = nlist.get(seed).nbr.get(j);
                // enqueue(), add to FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
                bfsQueue.enqueue(new int[]{seed, next});
                discovered[seed].set(j, true);                                 // set label to discovered in both neighbouting index lists
                discovered[next].set(nlist.get(next).nbr.indexOf(seed), true); // index where the background link was found
            }

            while (bfsQueue.hasItems()) {

                // dequeue(), take from FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
                int [] getLnk = (int[]) bfsQueue.dequeue();

                // next neighbour at the time it was added to the queue becomes current
                int prev = getLnk[0];
                int curr = getLnk[1];

                // always add the first node (it exists since this one was stored in the queue)
                Node n1 = new Node(nlist.get(curr).loc[0], nlist.get(curr).loc[1], nlist.get(curr).loc[2], nlist.get(curr).r, tree_count);
                n1.nbr.add(nodemap[prev]);
                nodemap[curr] = tree.size();
                tree.add(n1);

                while(Collections.frequency(discovered[curr], false)==1) { // while the number of undiscovered is 1 just step further

                    prev = curr;
                    curr = nlist.get(curr).nbr.get(discovered[curr].indexOf(false)); // curr takes the value of the next step

                    Node n2 = new Node(nlist.get(curr).loc[0], nlist.get(curr).loc[1], nlist.get(curr).loc[2], nlist.get(curr).r, tree_count);
                    n2.nbr.add(nodemap[prev]);
                    nodemap[curr] = tree.size();
                    tree.add(n2);

                    // mark as discovered the connections curr--prev and prev--curr
                    discovered[curr].set(nlist.get(curr).nbr.indexOf(prev), true);
                    discovered[prev].set(nlist.get(prev).nbr.indexOf(curr), true);

                }

                // means that it was not on the neurite anymore, check adjacent traces, now it is either endpoint or junction
                for (int i = 0; i < discovered[curr].size(); i++) {
                    boolean isDiscovered =  (Boolean) discovered[curr].get(i);
                    if (!isDiscovered) { // if it was not discovered
                        int next = nlist.get(curr).nbr.get(i);

                        bfsQueue.enqueue(new int[]{curr, next});   // enqueue()

                        discovered[curr].set(i, true); // label as discovered
                        discovered[next].set(nlist.get(next).nbr.indexOf(curr), true);

                    }
                }

            }

            tree_count++;

        }

        IJ.log(" -> "+tree_count+ " trees, "+ (tree.size()-1) +" nodes.");

        return tree;

    }

    class BfsQueue<E> {
        private LinkedList<E> list = new LinkedList<E>();
        public void enqueue(E item) {
            list.addLast(item);
        }
        public E dequeue() {
            return list.poll();
        }
        public boolean hasItems() {
            return !list.isEmpty();
        }
        public int size() {
            return list.size();
        }
//    public void addItems(BfsQueue<? extends E> q) {
//        while (q.hasItems())
//            list.addLast(q.dequeue());
//    }
    }

}
