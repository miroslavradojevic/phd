package com.braincadet.phd;

import com.braincadet.vess.Frangi;
//import features.TubenessProcessor;
import ij.*;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.plugin.filter.MaximumFinder;
import ij.process.*;

import java.awt.*;
import java.io.*;
import java.util.*;

public class MTracker implements PlugIn {

    String sigmas = "2,4";              // comma separated scale values (tubularity)

    // one sigmas tubularity image can be called for sequence of parameters (comma separated values)
    // parameter lists (string + array with extracted values)
    int[] no;                           // initial multi-object state cardinality
    String no_csv = "20";               // string with no values in CSV format

    int[] ro;                           // number of particles per object
    String ro_csv = "10";               // string with values in CSV format

    int[] ni;                           // number of predictions per particle
    String ni_csv = "10";               // string with ni

    int[] krad;                         // tube diameter in pixels
    String krad_csv = "4.0";

    int[] step;                         // motion step
    String step_csv = "3.0";

    float[] kappa;                      // von Mises angular probability
    String kappa_csv = "3";

    float[] pS;                         // survival pty
    String pS_csv = "0.95";

    float[] pD;                         // detection pty
    String pD_csv = "0.95";

    float[] th;                         // reference tubularity ratio for clutter
    String th_csv = "0.03";

    float[] kc;                         // clutter phd1 decay parameter
    String kc_csv = "4";

    int maxepoch = 150;                  // epoch limit
    int EPOCH_LOG = 10;                  // export reconstruction each EPOCH_LOG cycles

    // save results
    int maxiter = 200;                  // iteration limit (hundreds are fine)
    boolean savemidres = false;         // save partial results
    float[] tness;                      // tubeness min-max normalized

    int[] suppmap;                      // supression map: disable sampling (image stack size)
    int suppmap_width;
    int suppmap_height;
    int suppmap_length;

//    int N, M, P, SZ;                  // stack dimensions (width, height, length, size)
    String imdir, imnameshort;
    String midresdir = "";              // output directories, filenames

    // loggers (used in savemidres)
    public int[] X_cnt = new int[]{0}; // these are counters used in logging to count lines in swc output
    public String X_swclog = "";

    public int[] XP_cnt = new int[]{0}; //
    public String XP_swclog = "";

    public int[] ZP_cnt = new int[]{0}; // logging for the particles used to get the observations
    public String ZP_swclog = "";

    public int[] Z_cnt = new int[]{0};
    public String Z_swclog = "";

    public String zsizeCsvLog = "";       // iteration value at each particle
    public String tnessCsvLog = "";       // tness value per particle at each iteration
    public String phdmassCsvLog = "";       // phdmass at each particle

    int iter_count;                             // loop iterations

    private static final long MEGABYTE = 1024L * 1024L;

    private static long bytesToMegabytes(long bytes) {
        return bytes / MEGABYTE;
    }

    // simple struct containing image information
    private class Image8 {
        public int Width;
        public int Height;
        public int Length;
        String short_title;
        String image_dir;
        byte[] data; // byte8 data array
    }

    public void run(String s) {

        String image_path; // read input image, store the most recent path in Prefs
        String in_folder = Prefs.get("com.braincadet.phd.dir", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image");
        in_folder = dc.getDirectory();
        Prefs.set("com.braincadet.phd.dir", in_folder);
        image_path = dc.getPath();

        if (Macro.getOptions() == null) {
            GenericDialog gd = new GenericDialog("PHD");
            gd.addStringField("sigmas",     Prefs.get("com.braincadet.phd.sigmas", sigmas), 10);
            gd.addStringField("th",         Prefs.get("com.braincadet.phd.th", th_csv), 10);
            gd.addStringField("no",         Prefs.get("com.braincadet.phd.no", no_csv), 20);
            gd.addStringField("ro",         Prefs.get("com.braincadet.phd.ro", ro_csv), 10);
            gd.addStringField("ni",         Prefs.get("com.braincadet.phd.ni", ni_csv), 10);
            gd.addStringField("step",       Prefs.get("com.braincadet.phd.step", step_csv), 10);
            gd.addStringField("kappa",      Prefs.get("com.braincadet.phd.kappa", kappa_csv), 10);
            gd.addStringField("ps",         Prefs.get("com.braincadet.phd.ps", pS_csv), 10);
            gd.addStringField("pd",         Prefs.get("com.braincadet.phd.pd", pD_csv), 10);
            gd.addStringField("krad",       Prefs.get("com.braincadet.phd.krad", krad_csv), 10);
            gd.addStringField("kc",         Prefs.get("com.braincadet.phd.kc", kc_csv), 10);
            gd.addNumericField("maxiter",   Prefs.get("com.braincadet.phd.maxiter", maxiter), 0, 5, "");
            gd.addStringField("maxepoch",   Prefs.get("com.braincadet.phd.maxepoch", Integer.toString(maxepoch)), 10);
            gd.addCheckbox("savemidres",    Prefs.get("com.braincadet.phd.savemidres", savemidres));

            gd.showDialog();
            if (gd.wasCanceled()) return;

            sigmas = gd.getNextString();
            Prefs.set("com.braincadet.phd.sigmas", sigmas);
            th_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.th", th_csv);
            no_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.no", no_csv);
            ro_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.ro", ro_csv);
            ni_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.ni", ni_csv);
            step_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.step", step_csv);
            kappa_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.kappa", kappa_csv);
            pS_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.ps", pS_csv);
            pD_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.pd", pD_csv);
            krad_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.krad", krad_csv);
            kc_csv = gd.getNextString();
            Prefs.set("com.braincadet.phd.kc", kc_csv);
            maxiter = (int) gd.getNextNumber();
            Prefs.set("com.braincadet.phd.maxiter", maxiter);
            String maxepoch_str = gd.getNextString();
            maxepoch = (maxepoch_str.equals("inf")) ? Integer.MAX_VALUE : Integer.valueOf(maxepoch_str);
            Prefs.set("com.braincadet.phd.maxepoch", maxepoch);
            savemidres = gd.getNextBoolean();
            Prefs.set("com.braincadet.phd.savemidres", savemidres);
        } else {
            sigmas = Macro.getValue(Macro.getOptions(), "sigmas", sigmas);
            th_csv = Macro.getValue(Macro.getOptions(), "th", th_csv);
            no_csv = Macro.getValue(Macro.getOptions(), "no", no_csv);
            ro_csv = Macro.getValue(Macro.getOptions(), "ro", ro_csv);
            ni_csv = Macro.getValue(Macro.getOptions(), "ni", ni_csv);
            step_csv = Macro.getValue(Macro.getOptions(), "step", step_csv);
            kappa_csv = Macro.getValue(Macro.getOptions(), "kappa", kappa_csv);
            pS_csv = Macro.getValue(Macro.getOptions(), "ps", pS_csv);
            pD_csv = Macro.getValue(Macro.getOptions(), "pd", pD_csv);
            krad_csv = Macro.getValue(Macro.getOptions(), "krad", krad_csv);
            kc_csv = Macro.getValue(Macro.getOptions(), "kc", kc_csv);
            maxiter = Integer.valueOf(Macro.getValue(Macro.getOptions(), "maxiter", String.valueOf(maxiter)));
            String maxepoch_str = Macro.getValue(Macro.getOptions(), "maxepoch", String.valueOf(maxepoch));
            maxepoch = (maxepoch_str.equals("inf")) ? Integer.MAX_VALUE : Integer.valueOf(maxepoch_str);
            savemidres = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "savemidres", String.valueOf(savemidres)));
        }

//        image_path = "/Users/miroslav/exp.syn/test.65/65_cor=0.0_snr=5.tif";

        if (image_path == null) return;

        Image8 ip_load8 = load_image(image_path);

        if (ip_load8 == null) {
            IJ.log("could not load " + image_path);
            return;
        }

        int N = ip_load8.Width;  // ip_load.getWidth();
        int M = ip_load8.Height; // ip_load.getHeight();
        int P = ip_load8.Length; // ip_load.getStack().getSize();

        imnameshort = ip_load8.short_title; // ip_load.getShortTitle();
        imdir = ip_load8.image_dir; // ip_load.getOriginalFileInfo().directory;
        byte[] I = ip_load8.data;

//        // read image into byte[]
//        img = new float[SZ];
//        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1
//            byte[] slc = (byte[]) ip_load.getStack().getPixels(z);
//            for (int x = 0; x < N; x++) {
//                for (int y = 0; y < M; y++) {
//                    img[(z - 1) * (N * M) + y * N + x] = slc[y * N + x] & 0xff;
//                }
//            }
//        }

        if (savemidres) {
            midresdir =   ip_load8.image_dir + ip_load8.short_title + "_midres"; //  ip_load.getOriginalFileInfo().directory + ip_load.getTitle() + "_midres";
            Tools.createAndCleanDir(midresdir); // create midresult dir and initialize export/log
        }

        String[] dd;

        // extract comma separated parameter values
        dd = th_csv.split(",");
        if (dd.length == 0) return;
        th = new float[dd.length];
        for (int i = 0; i < dd.length; i++) th[i] = Float.valueOf(dd[i]);

        dd = no_csv.split(",");
        if (dd.length == 0) return;
        no = new int[dd.length];
        for (int i = 0; i < dd.length; i++) no[i] = Integer.valueOf(dd[i]);

        dd = ro_csv.split(",");
        if (dd.length == 0) return;
        ro = new int[dd.length];
        for (int i = 0; i < dd.length; i++) ro[i] = Integer.valueOf(dd[i]);

        dd = ni_csv.split(",");
        if (dd.length == 0) return;
        ni = new int[dd.length];
        for (int i = 0; i < dd.length; i++) ni[i] = Integer.valueOf(dd[i]);

        dd = step_csv.split(",");
        if (dd.length == 0) return;
        step = new int[dd.length];
        for (int i = 0; i < dd.length; i++) step[i] = Math.round(Float.valueOf(dd[i]));

        dd = kappa_csv.split(",");
        if (dd.length == 0) return;
        kappa = new float[dd.length];
        for (int i = 0; i < dd.length; i++) kappa[i] = Float.valueOf(dd[i]);

        dd = pS_csv.split(",");
        if (dd.length == 0) return;
        pS = new float[dd.length];
        for (int i = 0; i < dd.length; i++) pS[i] = Float.valueOf(dd[i]);

        dd = pD_csv.split(",");
        if (dd.length == 0) return;
        pD = new float[dd.length];
        for (int i = 0; i < dd.length; i++) pD[i] = Float.valueOf(dd[i]);

        dd = krad_csv.split(",");
        if (dd.length == 0) return;
        krad = new int[dd.length];
        for (int i = 0; i < dd.length; i++) krad[i] = Math.round(Float.valueOf(dd[i]));

        dd = kc_csv.split(",");
        if (dd.length == 0) return;
        kc = new float[dd.length];
        for (int i = 0; i < dd.length; i++) kc[i] = Float.valueOf(dd[i]);

        //******************************************************************
        // SYNTH: extract cor and snr parameters from the image name
//        dd = imnameshort.split("_");
//        IJ.log(Arrays.toString(dd));
//        imnameshort = dd[0];
//        float cor = Float.valueOf(dd[1].substring(4));
//        float snr = Float.valueOf(dd[2].substring(4));

        //******************************************************************
        IJ.log(" -- prefiltering...");

        // explanation ip_load is not necessary from this moment on, it was used for the fiji's tubularity
        // at this point ImagePlus needs to become byte[] ip_load_array, and ImagePlus object deleted..
        // since Frangi class uses byte[] as input and fiji's tubularity was using ImagePlus object

        long t1prep = System.currentTimeMillis();

//        is_tness = new ImageStack(N, M);
//        for (int i = 0; i < P; i++) {
//            float[] tt = new float[N * M];
//            Arrays.fill(tt, 0f);
//            is_tness.addSlice(new FloatProcessor(N, M, tt));
//        }

//        ip_tness = new ImagePlus("tness", is_tness);

        ArrayList<Float> sigs = new ArrayList<Float>();
        String[] readLn = sigmas.trim().split(",");

        for (int i = 0; i < readLn.length; i++) {
            sigs.add(Float.valueOf(readLn[i].trim()).floatValue());
//            TubenessProcessor tp = new TubenessProcessor(sig, false);
//            ImagePlus result = tp.generateImage(ip_load);
//            ImageCalculator ic = new ImageCalculator();

            // average, multipy
//          IJ.run(result, "Multiply...", "value=" + IJ.d2s(1f/readLn.length,3) + " stack");
//          ic.run("Add 32-bit stack", ip_tness, result); // result of the addition is placed in ip_tness

            // max
//            ic.run("Max 32-bit stack", ip_tness, result);
        }

        // plot sigs
//        for (int i = 0; i < sigs.size(); i++) {
//            IJ.log("sigs["+i+"] = " + sigs.get(i));
//        }
//        ip_tness.setCalibration(null);
//        if (savemidres) {
//            ImagePlus temp = ip_tness.duplicate();
//            IJ.run(temp, "8-bit", ""); // convert to 8 bit before saving
//            IJ.saveAs(temp, "Zip", midresdir + File.separator + "tness," + sigmas + ".zip");
//        }

        Frangi vess_filt = new Frangi(sigs); // instantiate new vesselness filtering class with the default parameters

        // float[] tness formation...[0,1] value range
        if (P>1){
            // 3d
//            t1 = System.currentTimeMillis();
//            J8 = f.run3d_byte(I, N, M, P, Vx, Vy, Vz, true); // threaded
            tness = vess_filt.run3d_float(I, N, M, P, null, null, null, true);
//            t2 = System.currentTimeMillis();
//            IJ.log("t = " + IJ.d2s((t2 - t1) / 1000f,2) + " [sec]");
//            new ImagePlus(imnameshort+"_Vess",  array2imagestack(J, N, M, P)).show();
        }
        else {
            // 2d
//            t1 = System.currentTimeMillis();
//            J8 = f.run2d_byte(I, N, M, P, Vx, Vy, Vz); // no threading
//            J = f.run2d_float(I, N, M, P, Vx, Vy, Vz);
            tness = vess_filt.run2d_float(I, N, M, P, null, null, null);
//            t2 = System.currentTimeMillis();
//            IJ.log("t = " + IJ.d2s((t2 - t1) / 1000f,2) + " [sec]");
//            new ImagePlus(imnameshort+"_Vess",  array2imagestack(J, N, M, P)).show();
        }

        // tubeness min-max normalize and store in an array for later and extract locations in a separate array
//        float tnessmin = Float.POSITIVE_INFINITY;
//        float tnessmax = Float.NEGATIVE_INFINITY;
//        tness = new float[SZ];

//        int[][] locationXYZ = new int[SZ][3]; // random sampling weighted with the normalized tubeness as importance function

//        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

//            float[] slc_float = null;
//            byte[] slc_byte = null;

//            if (usetness)
//                slc_float = (float[]) ip_tness.getStack().getPixels(z);
//            else
//                slc_byte = (byte[]) ip_load.getStack().getPixels(z);

//            for (int x = 0; x < N; x++) {
//                for (int y = 0; y < M; y++) {
//                    int ii = (z - 1) * (N * M) + y * N + x;

//                    if (usetness)
//                        tness[ii] = slc_float[y * N + x];
//                    else
//                        tness[ii] = (float) (slc_byte[y * N + x] & 0xff);

//                    locationXYZ[ii][0] = x;
//                    locationXYZ[ii][1] = y;
//                    locationXYZ[ii][2] = (z - 1);

//                    if (tness[ii] < tnessmin) tnessmin = tness[ii];
//                    if (tness[ii] > tnessmax) tnessmax = tness[ii];

//                }
//            }
//        }

//        for (int i = 0; i < SZ; i++) {
//            tness[i] = (tnessmax - tnessmin > Float.MIN_VALUE) ? ((tness[i] - tnessmin) / (tnessmax - tnessmin)) : 0;
//        }

        long t2prep = System.currentTimeMillis();
        IJ.log("t_prep = " + IJ.d2s((t2prep - t1prep) / 1000f,2) + " [sec]");

//        ImagePlus ip_vess   = new ImagePlus(imnameshort+"_Vess",  array2imagestack(tness, N, M, P));

        suppmap         = new int[N * M * P]; // suppression map with node tags
        suppmap_width   = ip_load8.Width;
        suppmap_height  = ip_load8.Height;
        suppmap_length  = ip_load8.Length;

        // go through comma separated parameter values
        for (int i01 = 0; i01 < no.length; i01++) {
            for (int i02 = 0; i02 < ro.length; i02++) {
                for (int i03 = 0; i03 < ni.length; i03++) {
                    for (int i04 = 0; i04 < krad.length; i04++) {
                        for (int i05 = 0; i05 < step.length; i05++) {
                            for (int i06 = 0; i06 < kappa.length; i06++) {
                                for (int i07 = 0; i07 < pS.length; i07++) {
                                    for (int i08 = 0; i08 < pD.length; i08++) {
                                        for (int i09 = 0; i09 < th.length; i09++) {
                                            for (int i10 = 0; i10 < kc.length; i10++) {

                                                if (savemidres) {

                                                    X_swclog = midresdir + File.separator + "Xk.swc";
                                                    X_cnt[0] = 0;
                                                    Tools.cleanfile(X_swclog);

                                                    XP_swclog = midresdir + File.separator + "XPk.swc";
                                                    XP_cnt[0] = 0;
                                                    Tools.cleanfile(XP_swclog);

                                                    ZP_swclog = midresdir + File.separator + "ZPk.swc";
                                                    ZP_cnt[0] = 0;
                                                    Tools.cleanfile(ZP_swclog);

                                                    Z_swclog = midresdir + File.separator + "Zk.swc";
                                                    Z_cnt[0] = 0;
                                                    Tools.cleanfile(Z_swclog);

                                                    tnessCsvLog = midresdir + File.separator + "tness.log";
                                                    Tools.cleanfile(tnessCsvLog);
                                                    zsizeCsvLog = midresdir + File.separator + "zsize.log";
                                                    Tools.cleanfile(zsizeCsvLog);
                                                    phdmassCsvLog = midresdir + File.separator + "phdmass.log";
                                                    Tools.cleanfile(phdmassCsvLog);
                                                }

                                                ArrayList<Integer> locs = new ArrayList<Integer>();     // list of candidate locations for seed points
                                                ArrayList<Float> locsw = new ArrayList<Float>();       // weights assigned to each location

                                                for (int z = ((P <= 3) ? 1 : 2); z <= ((P <= 3) ? P : P - 1); z++) { // layer count, zcoord is layer-1
                                                    Polygon maxx;
                                                    MaximumFinder mf = new MaximumFinder();

//                                                    maxx = mf.getMaxima(ip_tness.getStack().getProcessor(z), th[i09], false);

                                                    FloatProcessor fp = array2processor(tness, z-1, N, M);
                                                    maxx = mf.getMaxima(fp, th[i09], false); // get FloatProcessor from float[]


                                                    for (int i = 0; i < maxx.npoints; i++) {
                                                        int ii = (z - 1) * (N * M) + maxx.ypoints[i] * N + maxx.xpoints[i];
                                                        locs.add(ii);
                                                        locsw.add((float) Math.pow(tness[ii], MultiTT.weight_deg));
                                                    }
                                                }

                                                if (locs.size() == 0) {
                                                    IJ.log("0 seed candidate locations. no initiation for tolerance=" + IJ.d2s(th[i09], 2));
                                                    continue; // try another param configuration
                                                }

                                                float locs_count = locs.size();

                                                if (savemidres) {
                                                    // convert before exporting
                                                    ArrayList<int[]> tt = new ArrayList<int[]>(locs.size());

                                                    for (int i = 0; i < locs.size(); i++) {
                                                        int x = locs.get(i) % N;
                                                        int z = locs.get(i) / (N * M);
                                                        int y = locs.get(i) / N - z * M;
                                                        tt.add(new int[]{x, y, z});
                                                    }

                                                    exportlocsxyz(tt, 0.3f, Node.VIOLET, midresdir, "seedpool_tolerance=" + IJ.d2s(th[i09], 6));

                                                    tt.clear();
                                                }

                                                IJ.log("-- initialize...");
                                                MultiTT mtt; // multi-object tracker
                                                mtt = new MultiTT(P == 1, no[i01], ro[i02], ni[i03], krad[i04], step[i05], kappa[i06], pS[i07], pD[i08], th[i09], kc[i10]);

                                                if (savemidres) {
                                                    mtt.exporttemplates(midresdir);
                                                    Tools.createAndCleanDir(midresdir + File.separator + "mmodel");
                                                    mtt.mm.getModel(midresdir + File.separator + "mmodel");
                                                }

                                                Arrays.fill(suppmap, 0); // reset suppression map, it will fill up as rounds advance

                                                //******************************************************************
                                                // output directory name contains parameters
                                                String delindir = imdir + "PHD.sig.th.no.ro.ni.krad.stp.kapa.ps.pd.kc.e_" +
                                                        sigmas                  + "_" +
                                                        IJ.d2s(th[i09], 2)      + "_" +
                                                        IJ.d2s(no[i01], 0)      + "_" +
                                                        IJ.d2s(ro[i02], 0)      + "_" +
                                                        IJ.d2s(ni[i03], 0)      + "_" +
                                                        IJ.d2s(krad[i04], 0)    + "_" +
                                                        IJ.d2s(step[i05], 0)    + "_" +
                                                        IJ.d2s(kappa[i06], 1)   + "_" +
                                                        IJ.d2s(pS[i07], 2)      + "_" +
                                                        IJ.d2s(pD[i08], 2)      + "_" +
                                                        IJ.d2s(kc[i10], 1)      + "_";

                                                // SYNTH: redefine the output directory name to include cor and snr
//                                                delindir = imdir + "PHD.cor.snr.sig.th.no.ro.ni.krad.stp.kapa.ps.pd.kc.e_" +
//                                                        IJ.d2s(cor, 1)          + "_" + // to be compliant with dir naming in vaa3d experiments
//                                                        IJ.d2s(snr, 0)          + "_" + // expect integer snrs, to be compliant with dir naming in vaa3d experiments
//                                                        sigmas                  + "_" +
//                                                        IJ.d2s(th[i09], 2)      + "_" +
//                                                        IJ.d2s(no[i01], 0)      + "_" +
//                                                        IJ.d2s(ro[i02], 0)      + "_" +
//                                                        IJ.d2s(ni[i03], 0)      + "_" +
//                                                        IJ.d2s(krad[i04], 0)    + "_" +
//                                                        IJ.d2s(step[i05], 0)    + "_" +
//                                                        IJ.d2s(kappa[i06], 1)   + "_" +
//                                                        IJ.d2s(pS[i07], 2)      + "_" +
//                                                        IJ.d2s(pD[i08], 2)      + "_" +
//                                                        IJ.d2s(kc[i10], 1)      + "_";



                                                IJ.log("-- multi-object filtering...");
                                                long t1 = System.currentTimeMillis();

                                                int epochcnt = 0;

                                                while (locs.size() > 0 && epochcnt < maxepoch) {

                                                    epochcnt++;

                                                    iter_count = 0;

                                                    for (int i = locs.size() - 1; i >= 0; i--) { // exclude filtering from covered voxels
                                                        if (suppmap[locs.get(i)] >= mtt.suppmap_limit) {
                                                            locs.remove(i);
                                                            locsw.remove(i);
                                                        }
                                                    }

                                                    if (mtt.verbose)
                                                        IJ.log(IJ.d2s((locs.size() / locs_count) * 100f, 1) + "%");

                                                    if (locs.size() == 0) {

                                                        IJ.log("locs.size() == 0");

                                                        //-- export reconstruction
                                                        String outdir = delindir + IJ.d2s(epochcnt, 0);
                                                        Tools.createDir(outdir); // create if nonexistent
                                                        String SwcHeader = generateSwcHeader( // concatenate the params used into swc header format
                                                                imnameshort,
                                                                sigmas,
                                                                no[i01],
                                                                ro[i02],
                                                                krad[i04],
                                                                step[i05],
                                                                kappa[i06],
                                                                pS[i07],
                                                                pD[i08],
                                                                th[i09],
                                                                kc[i10],
                                                                maxiter,
                                                                maxepoch,
                                                                savemidres
                                                        );
                                                        export_reconstruction(mtt.Y, outdir+File.separator+imnameshort, SwcHeader);

                                                        break; // go out of while loop
                                                    }

                                                    ArrayList<int[]> N_o = new ArrayList<int[]>();

                                                    //********** initialization **********//
                                                    mtt._init(no[i01], step[i05], locs, locsw, N_o, N, M, P, tness, suppmap); // , template_ovrly img,
                                                    if (N_o.size() == 0) {

                                                        IJ.log("|N_o|=0");

                                                        // export reconstruction
                                                        String outdir = delindir + IJ.d2s(epochcnt, 0);
                                                        Tools.createDir(outdir); // create if nonexistent

                                                        String SwcHeader = generateSwcHeader( // concatenate the params used into swc header format
                                                                imnameshort,
                                                                sigmas,
                                                                no[i01],
                                                                ro[i02],
                                                                krad[i04],
                                                                step[i05],
                                                                kappa[i06],
                                                                pS[i07],
                                                                pD[i08],
                                                                th[i09],
                                                                kc[i10],
                                                                maxiter,
                                                                maxepoch,
                                                                savemidres
                                                        );

                                                        export_reconstruction(mtt.Y, outdir+File.separator+imnameshort, SwcHeader);

                                                        break; // out of while loop
                                                    }

                                                    if (savemidres) {

                                                        exportlocsxyz(N_o, 5f, Node.RED, midresdir, "seeds,e=" + IJ.d2s(epochcnt, 0));

                                                        exportXYZW(mtt.Xk, midresdir, "XWinit,epoch=" + IJ.d2s(epochcnt, 0), Node.MAGENTA);      // phd1 weights
                                                        exportXYZVxyz(mtt.Xk, midresdir, "XVinit,epoch=" + IJ.d2s(epochcnt, 0), Node.BLUE);     // directions of particles
                                                        exportNodes(mtt.Y, midresdir, "Yinit,epoch=" + IJ.d2s(epochcnt, 0)); // estimations

                                                        logval(tnessCsvLog, mtt.cluttertness);
                                                        logval(tnessCsvLog, mtt.Xk);
                                                        logval(zsizeCsvLog, no[i01]);                     // nr. observations mtt.Zk.size()
                                                        logval(phdmassCsvLog, mtt.phdmass);

                                                    }

                                                    if (mtt.Xk.size() > 0) {

                                                        while (iter_count < maxiter) {

                                                            boolean iterok;

                                                            iterok = mtt._iter1(N, M, P, tness, suppmap); // multi-object detection video demo: , template_ovrly

                                                            if (savemidres) { // iter_count%5==0 && iter_count==maxiter-1

                                                                Xlog(mtt.Xk, Node.GREEN);
                                                                XPlog(mtt.XPk, Node.BLUE_LIGHT);
                                                                Zlog(mtt.Zk, Node.RED, 1f);
                                                                ZPlog(mtt.ZPk, Node.OCRE, .1f);

                                                                logval(tnessCsvLog, mtt.XPk);
                                                                logval(zsizeCsvLog, mtt.Zk.size());
                                                                logval(phdmassCsvLog, mtt.phdmass);

                                                            }

                                                            if (!iterok)
                                                                break; // go out of while loop and try new set of seeds

                                                            iter_count++;

                                                        }

                                                    } else IJ.log("mtt.Xk.size() == 0");

                                                    if (locs.size() == 0 || epochcnt == maxepoch) { // put up with reconstruction every EPOCH_LOG: || epochcnt % EPOCH_LOG == 0
                                                        String outdir = delindir + IJ.d2s(epochcnt, 0);
                                                        Tools.createDir(outdir); // create if nonexistent

                                                        String SwcHeader = generateSwcHeader( // concatenate the params used into swc header format
                                                                imnameshort,
                                                                sigmas,
                                                                no[i01],
                                                                ro[i02],
                                                                krad[i04],
                                                                step[i05],
                                                                kappa[i06],
                                                                pS[i07],
                                                                pD[i08],
                                                                th[i09],
                                                                kc[i10],
                                                                maxiter,
                                                                maxepoch,
                                                                savemidres
                                                        );

                                                        export_reconstruction(mtt.Y, outdir+File.separator+imnameshort, SwcHeader);
                                                    }

                                                } // while there are locations and epochs have not reached the limit

                                                long t2 = System.currentTimeMillis();
                                                IJ.log("done. " + IJ.d2s(((t2 - t1) / 1000f), 2) + "s. [maxiter=" + maxiter + ", rounds=" + maxepoch + "]");

                                                // clear mtt components
                                                mtt.Xk.clear();
                                                mtt.Y.clear();

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

    private ImageStack array2imagestack(float[] a, int W, int H, int L) {

        ImageStack is = new ImageStack(W, H);

        for (int z = 0; z < L; z++) {
            float[] ai = new float[W*H];
            for (int j = 0; j < W * H; j++) {
                int x = j % W;
                int y = j / W;
                ai[j] = a[z*W*H+y*W+x];
            }

            is.addSlice(new FloatProcessor(W, H, ai));

        }

        return is;

    }

    private static String getFileExtension(String file_path)
    {
        String extension = "";

        int i = file_path.lastIndexOf('.');
        if (i >= 0) {
            extension = file_path.substring(i+1);
        }

        return extension;
    }

    private Image8 load_image(String image_path){

        if (
                !(
                getFileExtension(image_path).equalsIgnoreCase("tif") ||
                getFileExtension(image_path).equalsIgnoreCase("zip")
                )
                )
        {
            IJ.log("Input image can be tif or zip. Exiting...");
            return null;
        }

        ImagePlus input_image = new ImagePlus(image_path);

        if (input_image == null) {
            IJ.log("input_image == null");
            return null;
        }

        if (input_image.getType() != ImagePlus.GRAY8) {
            IJ.log("input_image.getType() != ImagePlus.GRAY8");
            return null;
        }

        IJ.run(input_image, "Flip Vertically", "stack"); // SYNTH: flip the tif image vertically before the computation, to be compliant with the vaa3d image matrix

        Image8 img = new Image8();
        img.Width   = input_image.getWidth();
        img.Height  = input_image.getHeight();
        img.Length  = input_image.getStack().getSize();
        img.short_title = input_image.getShortTitle();
        img.image_dir = input_image.getOriginalFileInfo().directory;

        img.data = new byte[img.Width*img.Height*img.Length]; // read image into byte[]
        for (int z = 1; z <= img.Length; z++) { // layer count, zcoord is layer-1
            byte[] slc = (byte[]) input_image.getStack().getPixels(z);
            for (int x = 0; x < img.Width; x++) {
                for (int y = 0; y < img.Height; y++) {
                    img.data[(z - 1) * (img.Width * img.Height) + y * img.Width + x] = slc[y * img.Width + x];// & 0xff;
                }
            }
        }

        return img;

    }

    private FloatProcessor array2processor(float[] a, int z_coord, int W, int H) {

        float[] lay = new float[W*H];

        for (int j = 0; j < W * H; j++) {
            int x = j % W;
            int y = j / W;
            lay[j] = a[z_coord *W*H + y*W + x];
        }

        return new FloatProcessor(W, H, lay);

    }

    private void export_reconstruction(
            ArrayList<Node> n0,
            String path_prefix,
            String SwcHeader
    )
    {

        float rad2kernel = 2f;
        int refine_iter = 3;
        double epsilon2 = 1e-8;
        float group_radius = 1.5f;

//        save_nodelist(n0, prefix+"_Phd0.swc", Node.RED);

        ArrayList<Node> n1 = new ArrayList<Node>();
        ArrayList<Node> n2 = new ArrayList<Node>();
        ArrayList<Node> n3 = new ArrayList<Node>();

        float resample_step = 1f;
        interpolate_nodelist(n0, resample_step, n1);
        non_blurring(n1, n2, rad2kernel, refine_iter, epsilon2);
        n1.clear();

//        int n_refinements = 1;
//        while (n_refinements < refine_iter) {
//            interpolate_nodelist(n2, resample_step, n1);
//            non_blurring(n1, n2, rad2kernel, 1, epsilon2);
//        }

        group(n2, n3, group_radius); // will do the connectivity fix
        n2.clear();

        ArrayList<Node> n3tree = bfs2(n3, true);
        n3.clear();

        save_nodelist(n3tree, path_prefix+".swc", Node.GREEN_LIGHT, SwcHeader);

    }

    private void save_nodelist(ArrayList<Node> nX, String SwcName, int type, String SwcHeader) {

        Tools.cleanfile(SwcName);

        try {

            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(SwcName, true)));

            logWriter.println(SwcHeader); // add swc header

            for (int i = 1; i < nX.size(); i++) { // nX[0] = null
                Node nn = nX.get(i);

                if (nn.nbr.size()==0) {
                    logWriter.println(IJ.d2s(i, 0)+" "+IJ.d2s((type<0)?nn.type:type, 0)+" "+IJ.d2s(nn.loc[0], 4)+" "+IJ.d2s(nn.loc[1], 4)+" "+IJ.d2s(nn.loc[2], 4)+" "+IJ.d2s(nn.r, 3)+" "+"-1");
                }
                else {
                    for (int j = 0; j < nn.nbr.size(); j++) {
                        logWriter.println(IJ.d2s(i, 0)+" "+IJ.d2s((type<0)?nn.type:type, 0)+" "+IJ.d2s(nn.loc[0], 4)+" "+IJ.d2s(nn.loc[1], 4)+" "+IJ.d2s(nn.loc[2], 4)+" "+IJ.d2s(nn.r, 3)+" "+IJ.d2s(nn.nbr.get(j), 0));
                    }
                }
            }

            logWriter.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void interpolate_nodelist(ArrayList<Node> nX, float resample_step, ArrayList<Node> nY) {

        // bidirectional connections between the nodes, interpolate inter-node line with resample_step size
        ArrayList<ArrayList<Boolean>> chk = new ArrayList<ArrayList<Boolean>>(nX.size()); // bookmark interpolated pairs
        nY.clear(); // nX copy, initialize output nY

        chk.add(null);
        nY.add(null);

        for (int i = 1; i < nX.size(); i++) {
            chk.add(new ArrayList<Boolean>(Collections.nCopies(nX.get(i).nbr.size(), false)));
            nY.add(new Node(nX.get(i)));
        }

        int init_size = nX.size();
        
        for (int i = 1; i < init_size; i++) {
            for (int j = 0; j < nX.get(i).nbr.size(); j++) {
                if (!chk.get(i).get(j)) {

                    int i1 = nX.get(i).nbr.get(j);
                    int j1 = nX.get(i1).nbr.indexOf(i);

                    if (j1 != -1) { // interpolate if there was existing link back

                        chk.get(i).set(j, true);
                        chk.get(i1).set(j1, true);

                        float vnorm = (float) Math.sqrt(
                                Math.pow(nX.get(i1).loc[0]-nX.get(i).loc[0], 2) +
                                Math.pow(nX.get(i1).loc[1]-nX.get(i).loc[1], 2) +
                                Math.pow(nX.get(i1).loc[2]-nX.get(i).loc[2], 2)
                        );

                        float vx = (nX.get(i1).loc[0]-nX.get(i).loc[0])/vnorm;
                        float vy = (nX.get(i1).loc[1]-nX.get(i).loc[1])/vnorm;
                        float vz = (nX.get(i1).loc[2]-nX.get(i).loc[2])/vnorm;
                        int N = (int) Math.round(vnorm/resample_step);

                        for (int k = 1; k < N; ++k) {
                            // add the node,only location is used in the refinement stage currently
                            nY.add(new Node(
                                    nX.get(i).loc[0] + k * (vnorm/N) * vx,
                                    nX.get(i).loc[1] + k * (vnorm/N) * vy,
                                    nX.get(i).loc[2] + k * (vnorm/N) * vz,
                                    nX.get(i).r + (nX.get(i1).r - nX.get(i).r) * (k/(float)N),
                                    ((k<=N/2)?nX.get(i).type:nX.get(i1).type)
                            ));

                            // backward
                            if (k==1) { // first
                                nY.get(nY.size()-1).nbr.add(i);
                                nY.get(i).nbr.set(j, nY.size()-1);
                            }
                            else { // middle
                                nY.get(nY.size()-1).nbr.add(nY.size()-2);
                                nY.get(nY.size()-2).nbr.add(nY.size()-1);
                            }

                            // forward
                            if (k==N-1) { // last
                                nY.get(nY.size()-1).nbr.add(i1);
                                nY.get(i1).nbr.set(j1, nY.size()-1);
                            }

                        }

                    }

                }
            }
        }

    }

    private void interpolate_treelist(ArrayList<Node> ntree, float resample_step, int type) {

        // assume onedirectional connections (1 link between 2 nodes in 1 direction)

        // interpolate all inter-node connections with the step size
        // not necessary to have bookkeeping variable, as there are no bidirectional links
        int init_size = ntree.size();
        for (int i = 1; i < init_size; i++) {
            // change the type
            if (type>=0) {
                // assign non-soma nodes with the given type
                if (ntree.get(i).type != Node.SOMA) {
                    ntree.get(i).type = type;
                }
            }

            // interpolation
            for (int j = 0; j < ntree.get(i).nbr.size(); ++j) { // there should be 0 or 1 neighbor
                int i1 = ntree.get(i).nbr.get(j);

                // interpolate between ntree[i] and ntree[i1]
                float vnorm = (float) Math.sqrt(Math.pow(ntree.get(i1).loc[0]-ntree.get(i).loc[0], 2) +
                        Math.pow(ntree.get(i1).loc[1]-ntree.get(i).loc[1], 2) +
                        Math.pow(ntree.get(i1).loc[2]-ntree.get(i).loc[2], 2));
                float vx = (ntree.get(i1).loc[0]-ntree.get(i).loc[0])/vnorm;
                float vy = (ntree.get(i1).loc[1]-ntree.get(i).loc[1])/vnorm;
                float vz = (ntree.get(i1).loc[2]-ntree.get(i).loc[2])/vnorm;
                int N = (int) Math.ceil(vnorm/resample_step);

                // add subsampling if N>1
                for (int k = 1; k < N; ++k) {
                    // add the node,only location is used in the refinement stage currently
                    ntree.add(new Node(
                            ntree.get(i).loc[0]+ k * (vnorm/N)*vx,
                            ntree.get(i).loc[1]+ k * (vnorm/N)*vy,
                            ntree.get(i).loc[2]+ k * (vnorm/N)*vz,
                            ntree.get(i).r  + (ntree.get(i1).r -ntree.get(i).r) * (k/(float)N),
                            ((k<=N/2)?ntree.get(i).type:ntree.get(i1).type)));

                    // link backward
                    if (k==1) {
                        // first
                        ntree.get(i).nbr.set(j, ntree.size()-1);
                    }
                    else {
                        // middle
                        ntree.get(ntree.size()-2).nbr.add(ntree.size()-1);
                    }
                    // link forward
                    if (k==N-1) {
                        // last
                        ntree.get(ntree.size()-1).nbr.add(i1);
                    }

                }

            }

        }

    }

    private void non_blurring(ArrayList<Node> nX, ArrayList<Node> nY, float rad2kernel, int MAXITER, double EPSILON2) {

        // mean-shift (non-blurring) uses flexible neighbourhood scaled with respect to the node's sigma
        int checkpoint = (int) Math.round(nX.size()/10.0);

        float[] conv = {0,0,0,0}; // x y z r
        float[] next = {0,0,0,0}; // x y z r

        nY.clear();
        // make a copy of the contents from nX, nY=nX
        for (int i = 0; i < nX.size(); i++) {
            nY.add(nX.get(i));
        }

        double x2, y2, z2, d2, r2;
        int iter, cnt;

        IJ.log("tree formation (experimental)...");

        for (int i = 1; i < nY.size(); i++) { // nY[0] is null

            if (i%checkpoint==0) IJ.log("" + ((i/checkpoint)*10) + "% " );

            conv[0] = nX.get(i).loc[0];
            conv[1] = nX.get(i).loc[1];
            conv[2] = nX.get(i).loc[2];
            conv[3] = nX.get(i).r;

            iter = 0;

            do { // iteratively change conv[]

                cnt = 0;

                next[0] = 0;
                next[1] = 0;
                next[2] = 0;
                next[3] = 0;

                r2 = Math.pow(rad2kernel * conv[3], 2);

                for (int j = 1; j < nX.size(); j++) {
                    x2 = Math.pow(nX.get(j).loc[0]-conv[0],2);
                    if (x2 <= r2) {
                        y2 = Math.pow(nX.get(j).loc[1]-conv[1],2);
                        if (x2 + y2 <= r2) {
                            z2 = Math.pow(nX.get(j).loc[2]-conv[2],2);
                            if (x2 + y2 + z2 <= r2) {

                                next[0] += nX.get(j).loc[0];
                                next[1] += nX.get(j).loc[1];
                                next[2] += nX.get(j).loc[2];
                                next[3] += nX.get(j).r;
                                cnt++;

                            }
                        }
                    }

                }

                next[0] /= cnt;
                next[1] /= cnt;
                next[2] /= cnt;
                next[3] /= cnt;

                d2 = Math.pow(next[0]-conv[0],2) + Math.pow(next[1]-conv[1],2) + Math.pow(next[2]-conv[2],2);

                conv[0] = next[0]; // for the next iteration
                conv[1] = next[1];
                conv[2] = next[2];
                conv[3] = next[3];

                iter++;

            }
            while(iter<MAXITER && d2>EPSILON2);

            nY.get(i).loc[0] = conv[0];
            nY.get(i).loc[1] = conv[1];
            nY.get(i).loc[2] = conv[2];
            nY.get(i).r      = conv[3];

        } // go through nY[i], initiate with nX[i] values and refine by mean-shift averaging

        IJ.log("done.");

    }

    private void group(ArrayList<Node> nX, ArrayList<Node> nY, float rad) { // sphere grouping

        ArrayList<Integer> X2Y = new ArrayList<Integer>(Collections.nCopies(nX.size(), -1)); // define mapping
        X2Y.set(0, 0); // the first (null) node corresponds to the first node in nY

        nY.clear();
        nY.add(null); // first node is the dummy node - null value in nY as it is in nX

        for (int i = 1; i < nX.size(); i++) { // add the rest of the nodes
            if (X2Y.get(i) != -1) continue; // skip if it was grouped already

            X2Y.set(i, nY.size());
            Node nYi = new Node(nX.get(i)); // initialize with the existing node
            float grp_size = 1;

            float r2 = rad * rad;
            double d2;
            for (int j = 1; j < nX.size(); j++) { // check the remainder of the nodes that was not grouped
                if (j!=i && X2Y.get(j)==-1) {
                    d2 = Math.pow(nX.get(j).loc[0]-nX.get(i).loc[0],2);
                    if (d2<=r2) {
                        d2 += Math.pow(nX.get(j).loc[1]-nX.get(i).loc[1],2); // [j].y
                        if (d2<=r2) {
                            d2 += Math.pow(nX.get(j).loc[2]-nX.get(i).loc[2],2);
                            if (d2<=r2) {

                                X2Y.set(j, nY.size()); // [j]=nY.size();

                                for (int k = 0; k < nX.get(j).nbr.size(); ++k)
                                    nYi.nbr.add(nX.get(j).nbr.get(k));  // append the neighbours of the group members

                                // update local average with x,y,z,r elements from nX[j]
                                grp_size++;
                                float a = (grp_size-1)/grp_size;
                                float b = (1f/grp_size);
                                nYi.loc[0]   = a * nYi.loc[0]     + b * nX.get(j).loc[0];
                                nYi.loc[1]   = a * nYi.loc[1]     + b * nX.get(j).loc[1];
                                nYi.loc[2]   = a * nYi.loc[2]     + b * nX.get(j).loc[2];
                                nYi.r        = a * nYi.r          + b * nX.get(j).r;

                            }
                        }
                    }

                }

            }

            nYi.type = Node.AXON;
            nY.add(nYi);

        }

        // refractor the indexing to nY
        for (int i = 1; i < nY.size(); i++) {
            for (int j = 0; j < nY.get(i).nbr.size(); ++j) {
                nY.get(i).nbr.set(j, X2Y.get(nY.get(i).nbr.get(j)));
            }
        }

        // remove double- and self- linkages after grouping
        connectivity_check(nY);

    }

    private void connectivity_check(ArrayList<Node> nX) {
        // - clean the bidirectional connections
        for (int i = 1; i < nX.size(); i++) {
            // remove double neighbourhoods from the neighbour list
            Set<Integer> set = new HashSet<Integer>();
            set.addAll(nX.get(i).nbr);
            nX.get(i).nbr.clear();
            nX.get(i).nbr.addAll(set);
            // remove self linkages (both solutions are fine if there was removal of the double neighbours)
            nX.get(i).nbr.remove((Integer) i);
//            nX.get(i).nbr.removeAll(Collections.singleton((Integer) i));
        }
        // ensure linkings are bidirectional, add if not
        for (int i = 1; i < nX.size(); i++) { // index 0 is null
            for (int j = 0; j < nX.get(i).nbr.size(); j++) {
                boolean fnd = false;
                for (int k = 0; k < nX.get(  nX.get(i).nbr.get(j)  ).nbr.size(); k++) {
                    if (nX.get(  nX.get(i).nbr.get(j)  ).nbr.get(k) == i) {
                        fnd = true;
                        break;
                    }
                }

                if (!fnd) {
                    nX.get( nX.get(i).nbr.get(j) ).nbr.add(i); // enforce link
                }
            }
        }


    }

    private int[][] sample(int nsamples, float[] csw, int[][] tosample) {

        int[][] out = new int[nsamples][tosample[0].length];

        float totalmass = csw[csw.length - 1];

        // use systematic resampling
        int i = 0;

        float u1 = (totalmass / (float) nsamples) * new Random().nextFloat();

        for (int j = 0; j < nsamples; j++) {

            float uj = u1 + j * (totalmass / (float) nsamples);

            while (uj > csw[i]) {
                i++;
            }

            for (int k = 0; k < tosample[i].length; k++) {
                out[j][k] = tosample[i][k];
            }

        }


        return out;
    }

    private void Xlog(ArrayList<X> x, int color) {
        loggerX(x, color, X_swclog, X_cnt);
    }

    private void XPlog(ArrayList<X> x, int color) {
        loggerX(x, color, XP_swclog, XP_cnt);
    }

    private void ZPlog(ArrayList<X> z, int color, float radius) {
        loggerZ(z, color, ZP_swclog, ZP_cnt, radius);
    }

    private void loggerX(ArrayList<X> x, int color, String swcfilepath, int[] swcindex) {

        try {

            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcfilepath, true)));

            for (int i = 0; i < x.size(); i++) {

                logWriter.println((++swcindex[0]) + " " + color + " " +
                        IJ.d2s(x.get(i).x, 3) + " " +
                        IJ.d2s(x.get(i).y, 3) + " " +
                        IJ.d2s(x.get(i).z, 3) + " " +
                        IJ.d2s(Math.pow((3 * x.get(i).w) / (4 * 3.14), 1f / 3), 3) + " " + (-1));

            }

            logWriter.close();

        } catch (IOException e) {}

    }

    private void Zlog(ArrayList<X> z, int color, float radius) {
        loggerZ(z, color, Z_swclog, Z_cnt, radius);
    }

    private void loggerZ(ArrayList<X> z, int color, String swcfilepath, int[] swcindex, float radius) {

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
                        IJ.d2s(z.get(i).x, 3) + " " +
                        IJ.d2s(z.get(i).y, 3) + " " +
                        IJ.d2s(z.get(i).z, 3) + " " +
                        IJ.d2s(radius, 3) + " " + (-1));//z.get(i).w

//                }


            }

            logWriter.close();

        } catch (IOException e) {
        }

    }

    private void logval(String filepath, ArrayList<X> x) {

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(filepath, true)));
            String s = "";
            for (int i = 0; i < x.size(); i++) s += ((i == 0) ? "" : ",") + IJ.d2s(x.get(i).tness, 4);
            logWriter.println(s);
            logWriter.close();
        } catch (IOException e) {
        }

    }

    private void logval(String filepath, float value) {

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(filepath, true)));
            logWriter.println(IJ.d2s(value, 4) + "");
            logWriter.close();
        } catch (IOException e) {
        }

    }

    private void exportlocsxyz(ArrayList<int[]> locsxyz, float radius, int swctype, String outdir, String swcname) {

        String outfile = outdir + File.separator + swcname + ".swc";

        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            for (int i = 0; i < locsxyz.size(); i++) {
                logWriter.println((i + 1) + " " +
                        swctype + " " +
                        IJ.d2s(locsxyz.get(i)[0], 4) + " " +
                        IJ.d2s(locsxyz.get(i)[1], 4) + " " +
                        IJ.d2s(locsxyz.get(i)[2], 4) + " " +
                        IJ.d2s(radius, 1) + " " + -1);
            }

            logWriter.close();

        } catch (IOException e) {
        }

    }

    private void exportXYZVxyz(ArrayList<X> Xobj, String outputdir, String swcname, int swctype) {

        String outfile = outputdir + File.separator + swcname + ".swc";

        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

//            logWriter.println("#Xk,Yk,Z,Vxyz");

            int swci = 0;
            for (int i = 0; i < Xobj.size(); i++) {
                logWriter.println((++swci) + " " + swctype + " " + IJ.d2s(Xobj.get(i).x, 3) + " " + IJ.d2s(Xobj.get(i).y, 3) + " " + IJ.d2s(Xobj.get(i).z, 2) + " " + 0.25 + " " + -1);
                int swcroot = swci;
                float dd = 5; // Xobj.get(i).sig*5f;
                logWriter.println((++swci) + " " + swctype + " " +
                        IJ.d2s(Xobj.get(i).x + dd * Xobj.get(i).vx, 3) + " " +
                        IJ.d2s(Xobj.get(i).y + dd * Xobj.get(i).vy, 3) + " " +
                        IJ.d2s(Xobj.get(i).z + dd * Xobj.get(i).vz, 3) + " " + 0.25 + " " + swcroot);
            }

            logWriter.close();

//            IJ.log("exported : " + outfile);

        } catch (IOException e) {
        }

    }

    private void exportXYZW(ArrayList<X> Xlist, String outdir, String swcname, int type) {

        String outfile = outdir + File.separator + swcname + ".swc";
        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

//            logWriter.println("#Xk,Yk,Z,W");

            for (int i = 0; i < Xlist.size(); i++) {
                logWriter.println((i + 1) + " " + type + " " +
                        IJ.d2s(Xlist.get(i).x, 4) + " " +
                        IJ.d2s(Xlist.get(i).y, 4) + " " +
                        IJ.d2s(Xlist.get(i).z, 4) + " " +
                        IJ.d2s(Math.pow((3 * Xlist.get(i).w) / (4 * 3.14), 1f / 3), 4) + " " + -1);
            }

            logWriter.close();

        } catch (IOException e) {
        }

    }

    private void exportNodes(ArrayList<Node> nlist, String outdir, String outname) {

        Tools.createDir(outdir);
        String recswc = outdir + File.separator + outname + ".swc";
        Tools.cleanfile(recswc);

        try {

            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(recswc, true)));

            int t = 0;

            for (int i = 0; i < nlist.size(); i++) {
                if (nlist.get(i) != null) {
                    logWriter.println((++t) + " " + ((i + 1) % 10) + " " + // YELLOW
                            IJ.d2s(nlist.get(i).loc[0], 3) + " " +
                            IJ.d2s(nlist.get(i).loc[1], 3) + " " +
                            IJ.d2s(nlist.get(i).loc[2], 3) + " " +
                            IJ.d2s(nlist.get(i).r, 3) + " " + (-1));
                }
            }

            logWriter.close();

        } catch (IOException e) {
        }

    }

    private void exportReconstruction(int numepochs, int numnodes, ArrayList<Node> nlist, String outdir, String outname) {//

        outdir += IJ.d2s(numepochs, 0) + "_" + IJ.d2s(numnodes, 0);
        Tools.createDir(outdir);
        String recswc = outdir + File.separator + outname + ".swc";
        Tools.cleanfile(recswc);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(recswc, true)));

            for (int i = 0; i < nlist.size(); i++) {
                if (nlist.get(i) != null) {

                    Node nn = nlist.get(i);

                    String out =
                            IJ.d2s(i, 0) + " " +
                                    IJ.d2s(nn.type, 0) + " " +
                                    IJ.d2s(nn.loc[0], 3) + " " +
                                    IJ.d2s(nn.loc[1], 3) + " " +
                                    IJ.d2s(nn.loc[2], 3) + " " +
                                    IJ.d2s(nn.r, 3) + " " +
                                    ((nn.nbr.size() == 0) ? "-1" : IJ.d2s(nn.nbr.get(0), 0));

                    logWriter.println(out);

                    if (nn.nbr.size() > 1)
                        IJ.log("*** ERROR in tree export " + i);
                }
            }

            logWriter.close();

        } catch (IOException e) {
        }

    }

    public ImagePlus getSuppMap() {

        ImageStack outis = new ImageStack(suppmap_width, suppmap_height);

        for (int z = 1; z <= suppmap_length; z++) { // layer count, zcoord is layer-1

            float[] slc = new float[suppmap_width * suppmap_height];

            for (int x = 0; x < suppmap_width; x++) {
                for (int y = 0; y < suppmap_height; y++) {

                    int ii = (z - 1) * (suppmap_width * suppmap_height) + y * suppmap_width + x;

                    slc[y * suppmap_width + x] = suppmap[ii];

                }
            }

            outis.addSlice(new FloatProcessor(suppmap_width, suppmap_height, slc));

        }

        return new ImagePlus("", outis); // title
    }

    private int[] gethist(ImagePlus imp, boolean usestack) {

        ImageProcessor ip = imp.getProcessor();

        ImageStatistics stats = null;
        if (usestack) {
            stats = new StackStatistics(imp);
        } else {
            if (!(ip instanceof ByteProcessor)) {
                ip.resetMinAndMax();
                imp.updateAndDraw();
            }
            stats = ImageStatistics.getStatistics(ip, Measurements.AREA + Measurements.MIN_MAX + Measurements.MODE, null);
        }

        int maxCount2 = 0;
        int[] histogram = new int[stats.nBins];

        for (int i = 0; i < stats.nBins; i++)
            histogram[i] = stats.histogram[i];

        for (int i = 0; i < stats.nBins; i++)
            if ((histogram[i] > maxCount2) && (i != stats.mode))
                maxCount2 = histogram[i];

        int hmax = stats.maxCount;

        if ((hmax > (maxCount2 * 2)) && (maxCount2 != 0)) {
            hmax = (int) (maxCount2 * 1.5);
            histogram[stats.mode] = hmax;
        }

        return histogram;

    }

    public void applythreshold(int thval, ImagePlus inimg) {

        int w = inimg.getWidth();
        int h = inimg.getHeight();
        int l = inimg.getStack().getSize();

        for (int i = 1; i <= l; i++) {
            byte[] lay = (byte[]) inimg.getStack().getProcessor(i).getPixels();
            for (int j = 0; j < w * h; j++) {
                if ((lay[j] & 0xff) >= thval) lay[j] = (byte) 255;
                else lay[j] = (byte) 0;
            }
        }

    }

    private String generateSwcHeader(
            String CurrentImname,
            String CurrentSigmas,
            int CurrentNo,
            int CurrentRo,
            int CurrentKrad,
            int CurrentStep,
            float CurrentKappa,
            float CurrentPs,
            float CurrentPd,
            float CurrentTh,
            float CurrentKc,
            int   CurrentMaxiter,
            int   CurrentMaxepoch,
            boolean CurrentSavemidres
    ){

        String outHeader = "#https://bitbucket.org/miroslavradojevic/phd/src\n";

        outHeader += "#imname = " + CurrentImname + "\n";
        outHeader += "#sigmas = " + CurrentSigmas + "\n";
        outHeader += "#no = "     + CurrentNo + "\n";
        outHeader += "#ro = "     + CurrentRo + "\n";
        outHeader += "#krad = "   + CurrentKrad + "\n";
        outHeader += "#step = "   + CurrentStep + "\n";
        outHeader += "#kappa = "  + CurrentKappa + "\n";
        outHeader += "#pS = "     + CurrentPs + "\n";
        outHeader += "#pD = "     + CurrentPd + "\n";
        outHeader += "#th = "     + CurrentTh + "\n";
        outHeader += "#Kc = "     + CurrentKc + "\n";
        outHeader += "#maxiter = " + CurrentMaxiter + "\n";
        outHeader += "#maxepoch = " + CurrentMaxepoch + "\n";
        outHeader += "#savemidres = " + CurrentSavemidres;

        return  outHeader;

    }

    private ArrayList<Integer> importsamp(ArrayList<Double> lcws, int n) {
        // systematic resampling, Beyond Kalman Filtering, Ristic et al.
        double totalmass = lcws.get(lcws.size() - 1);
        double u1 = (totalmass / (float) n) * new Random().nextDouble();

        ArrayList<Integer> out = new ArrayList<Integer>(n);
        out.clear();
        int i = 0;
        for (int j = 0; j < n; j++) {
            double uj = u1 + j * (totalmass / (float) n);
            while (uj > lcws.get(i)) i++;
            out.add(i);
        }

        return out;

    }

    private int[][] getoffsetsxyz(int rxy, int rz) { // call rz=0 for 2d image

        ArrayList<Integer> px = new ArrayList<Integer>();
        ArrayList<Integer> py = new ArrayList<Integer>();
        ArrayList<Integer> pz = new ArrayList<Integer>();

        for (int x = -rxy; x <= rxy; x++) {
            for (int y = -rxy; y <= rxy; y++) {
                for (int z = -rz; z <= rz; z++) {
                    if (Math.pow(x, 2) / Math.pow(rxy, 2) + Math.pow(y, 2) / Math.pow(rxy, 2) +
                            ((rz > 0) ? Math.pow(z, 2) / Math.pow(rz, 2) : 0)
                            <= 1) {
                        px.add(x);
                        py.add(y);
                        pz.add(z);
                    }
                }
            }
        }

//        IJ.log(px.size() + " points in the nbhood ellipsoid for Rxy="+rxy+", Rz="+rz);

        int[][] offxyz = new int[px.size()][3];

        for (int i = 0; i < px.size(); i++) {

            offxyz[i][0] = px.get(i);
            offxyz[i][1] = py.get(i);
            offxyz[i][2] = pz.get(i);

        }

        return offxyz;

    }

    private ArrayList<Node> bfs2(ArrayList<Node> nlist, boolean remove_isolated_node) {

        IJ.log("bfs...");

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
                if (nodesInTree==1 && !q.hasItems() && remove_isolated_node) {
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

    private class BfsQueue<E> {
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
    }

}
