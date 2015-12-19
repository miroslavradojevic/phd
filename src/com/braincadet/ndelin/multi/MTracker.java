package com.braincadet.ndelin.multi;

/*
Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

import com.braincadet.ndelin.fun.Tools;
import features.TubenessProcessor;
import ij.*;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class MTracker implements PlugIn {

    // color coding for the swc in vaa3d
    static int BLUE=3, YELLOW=6, RED=2, BLACK=1, MAGENTA=5, WHITE=0, VIOLET=4, GREEN=7, OCHRE=8, WEAK_GREEN=9, PINK=10, BLUEBERRY=12;

    static int NOBJ_INIT = 150;

    // allow having array of comma separated parameter inputs for the tubularity measure
    String  scales = "2,3";                   // comma separated scale values (tubularity)

    // one scales tubularity image can be called for sequence of parameters (comma separated values)
    // parameter lists (string + array with extracted values)
    int[]   nobjstart;                      // initial multi-object state cardinality
    String  nobjstart_csv = "";             // comma separated string of test parameters

    int[]   ro;                             // number of particles per object
    String  ro_csv = "";

    int[]   ni;                             // number of predictions per particle
    String  ni_csv = "";

    int[]   diam;                           // tube diameter in pixels
    String  diam_csv = "";

    int[]   step;                           // motion step
    String  step_csv = "";

    float[] kappa;                          // von Mises angular probability
    String  kappa_csv = "";

    float[] pS;                             // survival pty
    String  pS_csv = "";

    float[] pD;                             // detection pty
    String  pD_csv = "";

    float[] cluttertness;                   // reference tubularity ratio for clutter
    String  cluttertness_csv = "";

    float[] kclutt;                         // clutter phd decay parameter
    String  kclutt_csv = "";

    // save results
    int     maxiter = 100;                  // iteration limit
    int     logevery = 50;                  // iterations at which midresults are saved and the delineation output is stored
    boolean savemidres = false;             // save partial results

    float[] img;                            // original image as float array
    float[] tness;                          // tubeness min-max normalized
    int[]   suppmap;                          // supression map: disable sampling (image stack size)

//    int[][] locationXYZ;

    int N, M, P, SZ;                        // stack dimensions (width, height, length, size)
    String imdir, imnameshort;
    String midresdir;                       // output directories, filenames
    MultiTT mtt;                            // multi-object tracker

    // loggers
    public int[]    X_cnt       = new int[]{0};
    public String   X_swclog    = "";

    public int[]    XP_cnt      = new int[]{0};
    public String   XP_swclog   = "";

    public int[]    Z_cnt = new int[]{0};
    public String   Z_swclog = "";

    public String   zsizeCsvLog     = "";       // iteration value at each particle
    public String   tnessCsvLog     = "";       // tness value per particle at each iteration
    public String   phdmassCsvLog   = "";       // phdmass at each particle

    int iter_count;                             // loop iterations

    public void run(String s) {

        // read input image
        String in_folder = Prefs.get("com.braincadet.ndelin.multi.dir", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory(); Prefs.set("com.braincadet.ndelin.multi.dir", in_folder);
        String image_path = dc.getPath();
        if (image_path==null) return;

        ImagePlus ip_load = new ImagePlus(image_path);
        if(ip_load==null) {IJ.log(image_path + " was null"); return;}

        N = ip_load.getWidth();
        M = ip_load.getHeight();
        P = ip_load.getStack().getSize();
        SZ = N*M*P;

        ip_load.setCalibration(null);

        imnameshort = ip_load.getShortTitle();
        imdir = ip_load.getOriginalFileInfo().directory;

        // read image into byte[]
        img = new float[SZ];
        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1
            byte[] slc = (byte[]) ip_load.getStack().getPixels(z);
            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {
                    img[(z-1)*(N*M)+y*N+x] = slc[y*N+x]&0xff;
                }
            }
        }

        // load parameters
        if (Macro.getOptions()==null) {

            GenericDialog gd = new GenericDialog("TreeDelinPHD");

//            gd.addNumericField("nobjstart",             Prefs.get("com.braincadet.ndelin.multi.nobjstart", nobjstart), 0, 5, "["+NOBJ_INIT+"]");
//            gd.addNumericField("ro",                    Prefs.get("com.braincadet.ndelin.multi.ro", ro), 0, 5, "#");
//            gd.addNumericField("ni",                    Prefs.get("com.braincadet.ndelin.multi.ni", ni), 0, 5, "#");
//            gd.addNumericField("diam",                  Prefs.get("com.braincadet.ndelin.multi.diam", diam), 0, 5, "pix");
//            gd.addStringField("scales",                 Prefs.get("com.braincadet.ndelin.multi.scales", scales), 8);
//            gd.addNumericField("step",                  Prefs.get("com.braincadet.ndelin.multi.step", step), 0, 5, "pix");
//            gd.addNumericField("kappa",                 Prefs.get("com.braincadet.ndelin.multi.kappa", kappa), 1, 5, "(von Mises)");
//            gd.addNumericField("ps",                    Prefs.get("com.braincadet.ndelin.multi.ps", pS), 2, 5, "survival pty");
//            gd.addNumericField("pd",                    Prefs.get("com.braincadet.ndelin.multi.pd", pD), 2, 5, "detectionpty");
//            gd.addNumericField("cluttertness",          Prefs.get("com.braincadet.ndelin.multi.cluttertness", cluttertness), 2, 5, "[0-min tness, 1-begin tness]");
//            gd.addNumericField("kclutt",                Prefs.get("com.braincadet.ndelin.multi.kclutt", kclutt), 1, 5, "e^(-kclutt*tness)");

            gd.addStringField("scales",       Prefs.get("com.braincadet.ndelin.multi.scales", scales), 10);
            gd.addMessage("");
            gd.addStringField("nobjstart",    Prefs.get("com.braincadet.ndelin.multi.nobjstart", nobjstart_csv), 10);
            gd.addStringField("ro",           Prefs.get("com.braincadet.ndelin.multi.ro", ro_csv), 10);
            gd.addStringField("ni",           Prefs.get("com.braincadet.ndelin.multi.ni", ni_csv), 10);
            gd.addStringField("diam",         Prefs.get("com.braincadet.ndelin.multi.diam", diam_csv), 10);
            gd.addStringField("step",         Prefs.get("com.braincadet.ndelin.multi.step", step_csv), 10);
            gd.addStringField("kappa",        Prefs.get("com.braincadet.ndelin.multi.kappa", kappa_csv), 10);
            gd.addStringField("ps",           Prefs.get("com.braincadet.ndelin.multi.ps", pS_csv), 10);
            gd.addStringField("pd",           Prefs.get("com.braincadet.ndelin.multi.pd", pD_csv), 10);
            gd.addStringField("cluttertness", Prefs.get("com.braincadet.ndelin.multi.cluttertness", cluttertness_csv), 10);
            gd.addStringField("kclutt",       Prefs.get("com.braincadet.ndelin.multi.kclutt", kclutt_csv), 10);
            gd.addMessage("");
            gd.addNumericField("maxiter",     Prefs.get("com.braincadet.ndelin.multi.maxiter", maxiter), 0, 5, "");
            gd.addNumericField("logevery",    Prefs.get("com.braincadet.ndelin.multi.logevery", logevery), 0, 5, "");
            gd.addCheckbox("savemidres",      Prefs.get("com.braincadet.ndelin.multi.savemidres", savemidres));

            gd.showDialog();
            if (gd.wasCanceled()) return;

//            nobjstart   = (int) gd.getNextNumber();     Prefs.set("com.braincadet.ndelin.multi.nobjstart", nobjstart);
//            ro = (int) gd.getNextNumber();              Prefs.set("com.braincadet.ndelin.multi.ro", ro);
//            ni = (int) gd.getNextNumber();              Prefs.set("com.braincadet.ndelin.multi.ni", ni);
//            diam   = (int) gd.getNextNumber();          Prefs.set("com.braincadet.ndelin.multi.diam", diam);
//            step = (int) gd.getNextNumber();            Prefs.set("com.braincadet.ndelin.multi.step", step);
//            kappa = (float) gd.getNextNumber();         Prefs.set("com.braincadet.ndelin.multi.kappa", kappa);
//            pS = (float) gd.getNextNumber();            Prefs.set("com.braincadet.ndelin.multi.ps", pS);
//            pD = (float) gd.getNextNumber();            Prefs.set("com.braincadet.ndelin.multi.pd", pD);
//            cluttertness = (float) gd.getNextNumber();  Prefs.set("com.braincadet.ndelin.multi.cluttertness", cluttertness);
//            kclutt = (float) gd.getNextNumber();        Prefs.set("com.braincadet.ndelin.multi.kclutt", kclutt);

            scales      = gd.getNextString();     Prefs.set("com.braincadet.ndelin.multi.scales", scales);

            nobjstart_csv   = gd.getNextString(); Prefs.set("com.braincadet.ndelin.multi.nobjstart", nobjstart_csv);
            ro_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.ro", ro_csv);
            ni_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.ni", ni_csv);
            diam_csv   = gd.getNextString();      Prefs.set("com.braincadet.ndelin.multi.diam", diam_csv);
            step_csv = gd.getNextString();        Prefs.set("com.braincadet.ndelin.multi.step", step_csv);
            kappa_csv = gd.getNextString();       Prefs.set("com.braincadet.ndelin.multi.kappa", kappa_csv);
            pS_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.ps", pS_csv);
            pD_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.pd", pD_csv);
            cluttertness_csv = gd.getNextString();Prefs.set("com.braincadet.ndelin.multi.cluttertness", cluttertness_csv);
            kclutt_csv = gd.getNextString();      Prefs.set("com.braincadet.ndelin.multi.kclutt", kclutt_csv);

            maxiter = (int) gd.getNextNumber();   Prefs.set("com.braincadet.ndelin.multi.maxiter", maxiter);
            logevery = (int) gd.getNextNumber();  Prefs.set("com.braincadet.ndelin.multi.logevery", logevery);
            savemidres = gd.getNextBoolean();     Prefs.set("com.braincadet.ndelin.multi.savemidres", savemidres);

        }
        else {
            IJ.log("Macro.getOptions()");
//            scales = Macro.getValue(Macro.getOptions(),                    "scales", scales);IJ.log("scales="+scales);
//            nobjstart = Integer.valueOf(Macro.getValue(Macro.getOptions(), "nobjstart", String.valueOf(nobjstart)));IJ.log("nobjstart="+nobjstart);
//            ro = Integer.valueOf(Macro.getValue(Macro.getOptions(),        "ro", String.valueOf(ro)));IJ.log("ro="+ro);
//            ni = Integer.valueOf(Macro.getValue(Macro.getOptions(),        "ni", String.valueOf(ni))); IJ.log("ni="+ni);
//            diam = Integer.valueOf(Macro.getValue(Macro.getOptions(),      "diam", String.valueOf(diam)));IJ.log("diam="+diam);
//            step = Integer.valueOf(Macro.getValue(Macro.getOptions(),      "step", String.valueOf(step)));IJ.log("step="+step);
//            kappa = Float.valueOf(Macro.getValue(Macro.getOptions(),       "kappa", String.valueOf(kappa)));IJ.log("kappa="+kappa);
//            pS = Float.valueOf(Macro.getValue(Macro.getOptions(),            "ps", String.valueOf(pS)));IJ.log("ps="+pS);
//            pD = Float.valueOf(Macro.getValue(Macro.getOptions(),            "pd", String.valueOf(pD)));IJ.log("pd="+pD);
//            cluttertness = Float.valueOf(Macro.getValue(Macro.getOptions(),  "cluttertness", String.valueOf(cluttertness)));IJ.log("cluttertness="+cluttertness);
//            kclutt = Float.valueOf(Macro.getValue(Macro.getOptions(),        "kclutt", String.valueOf(kclutt))); IJ.log("kclutt="+kclutt);

            scales = Macro.getValue(Macro.getOptions(), "scales", scales);

            nobjstart_csv = Macro.getValue(Macro.getOptions(), "nobjstart", nobjstart_csv);
            ro_csv = Macro.getValue(Macro.getOptions(), "ro", ro_csv);
            ni_csv = Macro.getValue(Macro.getOptions(), "ni", ni_csv);
            diam_csv = Macro.getValue(Macro.getOptions(), "diam", diam_csv);
            step_csv = Macro.getValue(Macro.getOptions(), "step", step_csv);
            kappa_csv = Macro.getValue(Macro.getOptions(), "kappa", kappa_csv);
            pS_csv = Macro.getValue(Macro.getOptions(), "ps", pS_csv);
            pD_csv = Macro.getValue(Macro.getOptions(), "pd", pD_csv);
            cluttertness_csv = Macro.getValue(Macro.getOptions(), "cluttertness", cluttertness_csv);
            kclutt_csv = Macro.getValue(Macro.getOptions(), "kclutt", kclutt_csv);

            maxiter = Integer.valueOf(Macro.getValue(Macro.getOptions(),     "maxiter", String.valueOf(maxiter)));//IJ.log("maxiter="+maxiter);
            logevery = Integer.valueOf(Macro.getValue(Macro.getOptions(),    "logevery", String.valueOf(logevery)));//IJ.log("logevery="+logevery);
            savemidres = Boolean.valueOf(Macro.getValue(Macro.getOptions(),  "midres", String.valueOf(false)));//IJ.log("savemidres="+savemidres);
            IJ.log("");
        }

        String[] dd;

        // ge comma separated parameter values
        dd = nobjstart_csv.split(",");  if (dd.length==0) return;
        nobjstart = new int[dd.length];
        for (int i = 0; i < dd.length; i++) nobjstart[i] = Integer.valueOf(dd[i]);

        dd = ro_csv.split(","); if (dd.length==0) return;
        ro = new int[dd.length];
        for (int i = 0; i < dd.length; i++) ro[i] = Integer.valueOf(dd[i]);

        dd = ni_csv.split(","); if (dd.length==0) return;
        ni = new int[dd.length];
        for (int i = 0; i < dd.length; i++) ni[i] = Integer.valueOf(dd[i]);

        dd = diam_csv.split(","); if (dd.length==0) return;
        diam = new int[dd.length];
        for (int i = 0; i < dd.length; i++) diam[i] = Integer.valueOf(dd[i]);

        dd = step_csv.split(","); if (dd.length==0) return;
        step = new int[dd.length];
        for (int i = 0; i < dd.length; i++) step[i] = Integer.valueOf(dd[i]);

        dd = kappa_csv.split(","); if (dd.length==0) return;
        kappa = new float[dd.length];
        for (int i = 0; i < dd.length; i++) kappa[i] = Float.valueOf(dd[i]);

        dd = pS_csv.split(","); if (dd.length==0) return;
        pS = new float[dd.length];
        for (int i = 0; i < dd.length; i++) pS[i] = Float.valueOf(dd[i]);

        dd = pD_csv.split(","); if (dd.length==0) return;
        pD = new float[dd.length];
        for (int i = 0; i < dd.length; i++) pD[i] = Float.valueOf(dd[i]);

        dd = cluttertness_csv.split(","); if (dd.length==0) return;
        cluttertness = new float[dd.length];
        for (int i = 0; i < dd.length; i++) cluttertness[i] = Float.valueOf(dd[i]);

        dd = kclutt_csv.split(","); if (dd.length==0) return;
        kclutt = new float[dd.length];
        for (int i = 0; i < dd.length; i++) kclutt[i] = Float.valueOf(dd[i]);

        midresdir = ip_load.getOriginalFileInfo().directory+ip_load.getTitle()+"_midres"; // set output directory
        if (savemidres) {
            Tools.createAndCleanDir(midresdir); // create midresult dir and initialize export/log
            Tools.createAndCleanDir(midresdir + File.separator + "g(z|x)");
            Tools.createAndCleanDir(midresdir + File.separator + "suppmap");

            X_swclog  = midresdir    + File.separator + "Xk.swc";  X_cnt[0] = 0;    Tools.cleanfile(X_swclog);
            XP_swclog = midresdir    + File.separator + "XPk.swc"; XP_cnt[0] = 0;   Tools.cleanfile(XP_swclog);
            Z_swclog = midresdir    + File.separator +  "Zk.swc"; Z_cnt[0] = 0;     Tools.cleanfile(Z_swclog);

            tnessCsvLog = midresdir   + File.separator +    "tness.log";            Tools.cleanfile(tnessCsvLog);
            zsizeCsvLog = midresdir    + File.separator +   "zsize.log";            Tools.cleanfile(zsizeCsvLog);
            phdmassCsvLog = midresdir + File.separator +    "phdmass.log";          Tools.cleanfile(phdmassCsvLog);
        }

        //******************************************************************
        IJ.log(" -- tness");
        long t1 = System.currentTimeMillis();
        ImageStack  is_tness = new ImageStack(N, M);
        for (int i = 0; i < P; i++) {
            float[] tt = new float[N*M];
            is_tness.addSlice(new FloatProcessor(N, M, tt));
        }
        ImagePlus   ip_tness = new ImagePlus("tness", is_tness);
        String[] 	readLn = 	scales.trim().split(",");

        for (int i = 0; i < readLn.length; i++) {
            float sig = Float.valueOf(readLn[i].trim()).floatValue();
            TubenessProcessor tp = new TubenessProcessor(sig, false);
            ImagePlus result = tp.generateImage(ip_load);
            IJ.run(result, "Multiply...", "value=" + IJ.d2s(1f/readLn.length,3) + " stack");
            ImageCalculator ic = new ImageCalculator();
            ic.run("Add 32-bit stack", ip_tness, result); // result of the addition is placed in ip_tness
        }
        ip_tness.setCalibration(null);

        // tubeness min-max normalize and store in an array for later and extract locations in a separate array
        float tnessmin = Float.POSITIVE_INFINITY;
        float tnessmax = Float.NEGATIVE_INFINITY;
        tness = new float[SZ];
        int[][] locationXYZ = new int[SZ][3]; // random sampling weighted with the normalized tness as importance function

        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1
            float[] slc = (float[]) ip_tness.getStack().getPixels(z);
            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {
                    int ii = (z-1)*(N*M)+y*N+x;
                    tness[ii] = slc[y*N+x];
                    locationXYZ[ii][0] = x;
                    locationXYZ[ii][1] = y;
                    locationXYZ[ii][2] = (z-1);
                    if (tness[ii]<tnessmin) tnessmin = tness[ii];
                    if (tness[ii]>tnessmax) tnessmax = tness[ii];
                }
            }
        }
        long t2 = System.currentTimeMillis();
        IJ.log(((t2-t1)/1000f)+" sec.");


        if (savemidres) {
            IJ.run(ip_tness, "8-bit", ""); // convert to 8 bit before saving
            IJ.saveAs(ip_tness, "Tiff", midresdir + File.separator + "tness,"+scales+".tif");
        }

        //******************************************************************
        IJ.log(" -- sample tness");
        t1 = System.currentTimeMillis();
        int degree = 3;
        float[] p = new float[tness.length];
        for (int i = 0; i < tness.length; i++) {
            tness[i] = (tness[i]-tnessmin)/(tnessmax-tnessmin); // min-max normalize [0,1], will be used in importance sampling, directional filtering
            p[i] = (i==0)? (float) Math.pow(tness[0],degree) : p[i-1] + (float) Math.pow(tness[i], degree) ; // power to give more importance to those with high tness (todo: use cummulative weights)
        }
        int[][] samplexyz = sample(NOBJ_INIT, p, locationXYZ);
        t2 = System.currentTimeMillis();
        IJ.log(((t2-t1)/1000f)+" sec.");

        if (savemidres) {
            exportXYZR(samplexyz, 3.5f, VIOLET, midresdir, "Xsnap");
        }

        suppmap = new int[SZ];         // cummulative importance function (to avoid re-tracking)

        // go through comma separated parameter values
        for (int i01 = 0; i01 < nobjstart.length; i01++) {
            for (int i02 = 0; i02 < ro.length; i02++) {
                for (int i03 = 0; i03 < ni.length; i03++) {
                    for (int i04 = 0; i04 < diam.length; i04++) {
                        for (int i05 = 0; i05 <step.length; i05++) {
                            for (int i06 = 0; i06 < kappa.length; i06++) {
                                for (int i07 = 0; i07 < pS.length; i07++) {
                                    for (int i08 = 0; i08 < pD.length; i08++) {
                                        for (int i09 = 0; i09 < cluttertness.length; i09++) {
                                            for (int i10 = 0; i10 < kclutt.length; i10++) {

                                                String delindir = imdir+
                                                        "NDLN.scales.nobjstart.ro.ni.diam.step.kappa.ps.pd.cluttertness.kclutt.iter_"+
                                                        scales+"_"+
                                                        IJ.d2s(nobjstart[i01],0)+"_"+
                                                        IJ.d2s(ro[i02],0)+"_"+
                                                        IJ.d2s(ni[i03],0)+"_"+
                                                        IJ.d2s(diam[i04],0)+"_"+
                                                        IJ.d2s(step[i05],0)+"_"+
                                                        IJ.d2s(kappa[i06],1)+"_"+
                                                        IJ.d2s(pS[i07],2)+"_"+
                                                        IJ.d2s(pD[i08],2)+"_"+
                                                        IJ.d2s(cluttertness[i09],2)+"_"+
                                                        IJ.d2s(kclutt[i10],1)+"_";//+IJ.d2s(iteration,0);

                                                IJ.log(
                                                        "\nnobjstart="+nobjstart[i01]+" "+ "ro="+ro[i02]+" "+ "ni="+ni[i03]+" "+
                                                                "diam="+diam[i04]+"\n"+
                                                                "step="+step[i05]+"\n"+
                                                                "kappa="+kappa[i06]+"\n"+
                                                                "ps="+pS[i07]+" "+"pd="+pD[i08]+"\n"+
                                                                "cluttertness="+cluttertness[i09]+"\n"+
                                                                "kclutt="+kclutt[i10]+"\n"
                                                );

                                                //******************************************************************

                                                if (true) {

                                                    IJ.log("-- initialize...");
                                                    mtt = new MultiTT(
                                                            P == 1,
                                                            nobjstart[i01],
                                                            ro[i02],
                                                            ni[i03],
                                                            diam[i04],
                                                            step[i05],
                                                            kappa[i06],
                                                            pS[i07],
                                                            pD[i08],
                                                            cluttertness[i09],
                                                            kclutt[i10]);

                                                    if (savemidres) {
                                                        mtt.exporttemplates(midresdir);
                                                        Tools.createAndCleanDir(midresdir + File.separator + "mmodel");
                                                        mtt.mm.getModel(midresdir + File.separator + "mmodel");
                                                    }

                                                    //******************************************************************
                                                    Arrays.fill(suppmap, -2); // reset suppression map

                                                    //******************************************************************
                                                    IJ.log("-- multi-object filtering...");
                                                    t1 = System.currentTimeMillis();
                                                    iter_count = 0;
                                                    mtt.init(img, N, M, P, samplexyz, tness, suppmap);

                                                    //******************************************************************
                                                    if (savemidres) {
                                                        exportXYZW(mtt.Xk, midresdir, "X0_W", RED);      // phd weights
                                                        exportXYZVxyz(mtt.Xk, midresdir, "X0_V", BLUE);     // directions of particles
                                                        exportXYZSig(mtt.Xk, midresdir, "X0_Sig", MAGENTA);  // particles with their sigma
                                                        exportXYZSig(mtt.Yk, midresdir, "Y0", YELLOW);  // estimates

                                                        logval(tnessCsvLog, (float) mtt.tnessinit);
                                                        logval(tnessCsvLog, (float) mtt.tnessinit * mtt.cluttertness);
                                                        logval(tnessCsvLog, mtt.Xk);
                                                        logval(zsizeCsvLog, nobjstart[i01]);                     // nr. observations mtt.Zk.size()
                                                        logval(phdmassCsvLog, mtt.phdmass);                     // total phd mass log

                                                        String name = "suppmap,iter=" + IJ.d2s(-1, 0);
                                                        ImagePlus hpimp = getSuppMap(name);
                                                        IJ.saveAs(hpimp, "Tiff", midresdir + File.separator + "suppmap" + File.separator + name + ".tif");
                                                    }

                                                    while (iter_count < maxiter) { // && iter_count < mtt.ITER_LIMIT

                                                        boolean iterok = mtt.iter(iter_count, iter_count == 0, tness, N, M, P, tness, suppmap);

                                                        if (!iterok) {
                                                            break;
                                                        }

                                                        if (savemidres) {

                                                            Xlog(mtt.Xk, MAGENTA);
                                                            XPlog(mtt.XPk, WEAK_GREEN);
                                                            Zlog(mtt.Zk, RED);

                                                            String name = "suppmap,iter=" + IJ.d2s(iter_count, 0);
                                                            ImagePlus hpimp = getSuppMap(name);
                                                            IJ.saveAs(hpimp, "Tiff", midresdir + File.separator + "suppmap" + File.separator + name + ".tif");

                                                            logval(tnessCsvLog, mtt.XPk);
                                                            logval(zsizeCsvLog, mtt.Zk.size());
                                                            logval(phdmassCsvLog, mtt.phdmass);

                                                            ImagePlus gimp = new ImagePlus("g(z|x),iter=" + IJ.d2s(iter_count, 0), new FloatProcessor(mtt.g));
                                                            IJ.run(gimp, "Rotate 90 Degrees Right", "");
                                                            IJ.saveAs(gimp, "Tiff", midresdir + File.separator + "g(z|x)" + File.separator + "g(z|x),iter=" + IJ.d2s(iter_count, 0) + ".tif");
                                                        }

                                                        iter_count++;
                                                        if (iter_count % logevery == 0)
                                                            Yklog(iter_count, delindir, YELLOW, mtt.Y);

                                                    }

                                                    t2 = System.currentTimeMillis();
                                                    IJ.log("done. " + IJ.d2s((t2 - t1) / 1000f, 2) + "s." + " [maxiter=" + maxiter + "]");
                                                    Yklog(iter_count, delindir, YELLOW, mtt.Y);

                                                }

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // automatize thresholding in the code
//            // tness threshold global
//            IJ.run(ip_tness, "Smooth", "stack");
//            Object[] result = new Thresholder().exec(ip_tness, "Default", false, false, true, false, false, true);
//            ImagePlus iexapmle = (ImagePlus)result[1];
//            iexapmle.setTitle("observationmap");
//            if (savemidres) IJ.saveAs(iexapmle, "Tiff", midresdir + File.separator + "omap.tif");

    }

    private int[][] sample(int nsamples, float[] csw, int[][] tosample) {

        int[][] out = new int[nsamples][tosample[0].length];

        float totalmass = csw[csw.length-1];

        // use systematic resampling
        int i = 0;

        float u1 = (totalmass/(float)nsamples) * new Random().nextFloat();

        for (int j = 0; j < nsamples; j++) {

            float uj = u1 + j*(totalmass/(float)nsamples);

            while (uj > csw[i]) {
                i++;
            }

            for (int k = 0; k < tosample[i].length; k++) {
                out[j][k] = tosample[i][k];
            }

        }


        return out;
    }

    private void Yklog(int iteration, String delindir, int color, ArrayList<ArrayList<X>> Yk) {


        delindir += IJ.d2s(iteration,0);
        Tools.createDir(delindir);
        String delinswc = delindir + File.separator + imnameshort + ".swc";
        Tools.cleanfile(delinswc);

        try {

            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(delinswc, true)));

            int dd = 0;

            for (int i = 0; i < Yk.size(); i++) {

                for (int j = 0; j < Yk.get(i).size(); j++) {

                    X pcl = Yk.get(i).get(j);

                    logWriter.println((++dd) + " " + color + " " +
                            IJ.d2s(pcl.x,3) + " " +
                            IJ.d2s(pcl.y,3) + " " +
                            IJ.d2s(pcl.z,3) + " " +
                            IJ.d2s(1f,3) + " " + (-1));

                }

            }

            logWriter.close();

        } catch (IOException e) {}

        IJ.log(delinswc);

    }

    private void Xlog(ArrayList<X> x, int color) {
        loggerX(x, color, X_swclog, X_cnt);
    }

    private void XPlog(ArrayList<X> x, int color) {
        loggerX(x, color, XP_swclog, XP_cnt);
    }

    private void loggerX(ArrayList<X> x, int color, String swcfilepath, int[] swcindex) {
        
        try {

            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcfilepath, true)));

            for (int i = 0; i < x.size(); i++) {

                logWriter.println((++swcindex[0]) + " " + color + " " +
                        IJ.d2s(x.get(i).x,3) + " " +
                        IJ.d2s(x.get(i).y,3) + " " +
                        IJ.d2s(x.get(i).z,3) + " " +
                        IJ.d2s(Math.pow((3*x.get(i).w)/(4*3.14),1f/3),3) + " " + (-1));

            }

            logWriter.close();

        } catch (IOException e) {}
        
    }

    private void Zlog(ArrayList<X> z, int color) {
        loggerZ(z, color, Z_swclog, Z_cnt);
    }

    private void loggerZ(ArrayList<X> z, int color, String swcfilepath, int[] swcindex) {

        try {

            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcfilepath, true)));

            for (int i = 0; i < z.size(); i++) {
//                if (z.get(i).count==1) {
//
//                    logWriter.println((++swcindex[0]) + " " + BLUE + " " +
//                            IJ.d2s(z.get(i).x,3) + " " +
//                            IJ.d2s(z.get(i).y,3) + " " +
//                            IJ.d2s(z.get(i).z,3) + " " +
//                            IJ.d2s(z.get(i).sig,3) + " " + (-1));
//
//                }
//                else if (z.get(i).count==2) {
//
//                    logWriter.println((++swcindex[0]) + " " + BLUEBERRY + " " +
//                            IJ.d2s(z.get(i).x,3) + " " +
//                            IJ.d2s(z.get(i).y,3) + " " +
//                            IJ.d2s(z.get(i).z,3) + " " +
//                            IJ.d2s(z.get(i).sig,3) + " " + (-1));
//
//                }
//                else {

                    logWriter.println((++swcindex[0]) + " " + color + " " +
                            IJ.d2s(z.get(i).x,3) + " " +
                            IJ.d2s(z.get(i).y,3) + " " +
                            IJ.d2s(z.get(i).z,3) + " " +
                            IJ.d2s(1f,3) + " " + (-1));//z.get(i).w

//                }



            }

            logWriter.close();

        } catch (IOException e) {}

    }

    private void logval(String filepath, ArrayList<X> x) {

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(filepath, true)));
            String s = "";
            for (int i = 0; i < x.size(); i++) s += ((i == 0) ? "" : ",") + IJ.d2s(x.get(i).tness, 4);
            logWriter.println(s);
            logWriter.close();
        } catch (IOException e) {}

    }

    private void logval(String filepath, float value) {

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(filepath, true)));
            logWriter.println(IJ.d2s(value, 4) + "");
            logWriter.close();
        } catch (IOException e) {}

    }

    private void exportXYZR(int[][] locsxyz, float radius, int swctype, String outdir, String swcname) {

        String outfile = outdir + File.separator + swcname + ".swc";

        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            for (int i = 0; i < locsxyz.length; i++) {
                logWriter.println((i+1) + " " +
                        swctype + " " +
                        IJ.d2s(locsxyz[i][0], 4) + " " +
                        IJ.d2s(locsxyz[i][1], 4) + " " +
                        IJ.d2s(locsxyz[i][2], 4) + " " +
                        radius + " " + -1);
            }

            logWriter.close();

        } catch (IOException e) {}

    }

    private void exportXYZVxyz(ArrayList<X> Xobj, String outputdir, String swcname, int swctype) {

        String outfile = outputdir + File.separator + swcname + ".swc";

        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

//            logWriter.println("#Xk,Yk,Z,Vxyz");

            int swci = 0;
            for (int i = 0; i < Xobj.size(); i++) {
                logWriter.println((++swci) + " " + swctype + " " + IJ.d2s(Xobj.get(i).x,3) + " " + IJ.d2s(Xobj.get(i).y,3) + " " + IJ.d2s(Xobj.get(i).z,2) + " " + 0.25 + " " + -1);
                int swcroot = swci;
                float dd = Xobj.get(i).sig*5f;
                logWriter.println((++swci) + " " + swctype + " " + IJ.d2s(Xobj.get(i).x+ dd *Xobj.get(i).vx,3) + " " + IJ.d2s(Xobj.get(i).y+ dd *Xobj.get(i).vy,3) + " " + IJ.d2s(Xobj.get(i).z+ dd *Xobj.get(i).vz,3) + " " + 0.25 + " " + swcroot);
            }

            logWriter.close();

//            IJ.log("exported : " + outfile);

        } catch (IOException e) {}

    }

    private void exportXYZSig(ArrayList<X> Xlist, String outdir, String swcname, int type) {

        // exports only locations into swc format for visualization - the rest of the Xk instance is ignored
        if (outdir==null || swcname==null) return;

        String outfile = outdir + File.separator + swcname + ".swc";

        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            logWriter.println("#Xk,Yk,Z,Sig");

            for (int i = 0; i < Xlist.size(); i++) {
                logWriter.println((i+1) + " " + type + " " +
                        IJ.d2s(Xlist.get(i).x, 4) + " " +
                        IJ.d2s(Xlist.get(i).y, 4) + " " +
                        IJ.d2s(Xlist.get(i).z, 4) + " " +
                        Xlist.get(i).sig + " " + -1);
            }

            logWriter.close();

        } catch (IOException e) {}

    }

    private void exportXYZW(ArrayList<X> Xlist, String outdir, String swcname, int type) {

        String outfile = outdir + File.separator + swcname + ".swc";
        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            logWriter.println("#Xk,Yk,Z,W");

            for (int i = 0; i < Xlist.size(); i++) {
                logWriter.println((i+1) + " " + type + " " +
                        IJ.d2s(Xlist.get(i).x, 4) + " " +
                        IJ.d2s(Xlist.get(i).y, 4) + " " +
                        IJ.d2s(Xlist.get(i).z, 4) + " " +
                        IJ.d2s(Math.pow((3*Xlist.get(i).w)/(4*3.14),1f/3), 4) + " " + -1);
            }

            logWriter.close();

        } catch (IOException e) {}

    }

    public ImagePlus getSuppMap(String title){

        ImageStack outis = new ImageStack(N, M);

        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

            float[] slc = new float[N*M]; // (float[]) ip_tness.getStack().getPixels(z);

            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {

                    int ii = (z-1)*(N*M)+y*N+x;

                    slc[y*N+x] = suppmap[ii];

                }
            }

            outis.addSlice(new FloatProcessor(N, M, slc));

        }

        return new ImagePlus(title, outis);
    }

}
