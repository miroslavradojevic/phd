package com.braincadet.phd;

import com.braincadet.phd.fun.Tools;
import com.braincadet.phd.multi.MultiTT;
import com.braincadet.phd.multi.X;
import com.braincadet.phd.swc.Node;
import features.TubenessProcessor;
import ij.*;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.measure.Measurements;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.plugin.filter.MaximumFinder;
import ij.process.*;

import java.awt.*;
import java.io.*;
import java.util.*;

public class MTracker implements PlugIn {

    // color coding for the swc rendering as in Vaa3D swc visualization (www.vaa3d.org)
    static int WHITE = 0, BLACK = 1, RED = 2, BLUE = 3, VIOLET = 4, MAGENTA = 5, YELLOW = 6, GREEN = 7, OCHRE = 8, WEAK_GREEN = 9, PINK = 10, BLUEBERRY = 12;

    AutoThresholder.Method thmethod = AutoThresholder.Method.IsoData; // IsoData Moments Otsu Triangle Default

    // allow having array of comma separated parameter inputs for the tubularity measure
    String sigmas = "2,3";              // comma separated scale values (tubularity)

    // one sigmas tubularity image can be called for sequence of parameters (comma separated values)
    // parameter lists (string + array with extracted values)
    int[] no;              // initial multi-object state cardinality
    String no_csv = "20";   // string with no values in CSV format

    int[] ro;                 // number of particles per object
    String ro_csv = "10";      // string with values in CSV format

    int[] ni;                 // number of predictions per particle
    String ni_csv = "10";      // string with ni

    int[] krad;                         // tube diameter in pixels
    String krad_csv = "4.0";            //

    int[] step;                       // motion step
    String step_csv = "3.0";           //

    float[] kappa;                      // von Mises angular probability
    String kappa_csv = "3";

    float[] pS;                         // survival pty
    String pS_csv = "0.95";

    float[] pD;                         // detection pty
    String pD_csv = "0.95";

    float[] th;                         // reference tubularity ratio for clutter
    String th_csv = "10";

    float[] kc;                         // clutter phd decay parameter
    String kc_csv = "4";

    int maxepoch = 10;                  // epoch limit
//    String maxepoch_csv = "";         // comma separated input

    int TUBE_RADIUS = 3;                // used at _init(), radius of the sphere used for the seed point

    // save results
    int maxiter = 200;              // iteration limit (hundreds are fine)
    boolean savemidres = false;     // save partial results
    boolean usetness = true;        // use tubularity measure

    float[] img;                    // original image as float array
    float[] tness;                  // tubeness min-max normalized
    int[] suppmap;                  // supression map: disable sampling (image stack size)

    // template_* variables were used for the multi-object detection video demo
//    ImageStack template_stack;
//    ImagePlus template_image;
//    Overlay template_ovrly;
//    ZProjector template_zprojector;

    int N, M, P, SZ;                    // stack dimensions (width, height, length, size)
    String imdir, imnameshort;
    String midresdir = "";              // output directories, filenames

    // loggers
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

    public static long bytesToMegabytes(long bytes) {
        return bytes / MEGABYTE;
    }

    public void run(String s) {

        // read input image, store the most recent path in Prefs
        String in_folder = Prefs.get("com.braincadet.phd.dir", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image");
        in_folder = dc.getDirectory();
        Prefs.set("com.braincadet.phd.dir", in_folder);
        String image_path = dc.getPath();
        if (image_path == null) return;

        ImagePlus ip_load = new ImagePlus(image_path);

        if (ip_load == null) {
            IJ.log(image_path + " was null");
            return;
        }
        if (ip_load.getType() != ImagePlus.GRAY8) {
            IJ.log("Image needs to be GRAY8.");
            return;
        }

        N = ip_load.getWidth();
        M = ip_load.getHeight();
        P = ip_load.getStack().getSize();
        SZ = N * M * P;

        ip_load.setCalibration(null);

        imnameshort = ip_load.getShortTitle();
        imdir = ip_load.getOriginalFileInfo().directory;

//        if (false) {
//            // experimental code, prototype video export
//            Overlay o = new Overlay();
//            Color cc = new Color(1f,1f,0f,0.4f);
//            float rr = 1f;
//            for (int z = 1; z < P; z++) {
//                for (int x = 0; x < N/10; x++) {
//                    for (int y = 0; y < M/10; y++) {
//                        OvalRoi oi = new OvalRoi(x-rr+0.5, y-rr+0.5, 2*rr, 2*rr);
//                        oi.setPosition(z);
//                        oi.setFillColor(cc);
//                        oi.setStrokeColor(cc);
//                        o.add(oi);
//                    }
//                }
//            }
//
//            ImagePlus ip1 = ip_load.duplicate();
//            ip1.setOverlay(o);
//            ip1.updateAndDraw();
//            IJ.log("flatten...");
//            ip1.flattenStack();
//            FileSaver fs = new FileSaver(ip1);
//            fs.saveAsTiff(imdir+File.separator+imnameshort+"_flatt.tif");
//            IJ.log("saved!");
//            return;
//        }

        // read image into byte[]
        img = new float[SZ];
        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1
            byte[] slc = (byte[]) ip_load.getStack().getPixels(z);
            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {
                    img[(z - 1) * (N * M) + y * N + x] = slc[y * N + x] & 0xff;
                }
            }
        }

        if (Macro.getOptions() == null) {

            GenericDialog gd = new GenericDialog("PHD");
            gd.addStringField("sigmas", Prefs.get("com.braincadet.phd.sigmas", sigmas), 10);
            gd.addStringField("th", Prefs.get("com.braincadet.phd.th", th_csv), 10);
            gd.addStringField("no", Prefs.get("com.braincadet.phd.no", no_csv), 20);
            gd.addStringField("ro", Prefs.get("com.braincadet.phd.ro", ro_csv), 10);
            gd.addStringField("ni", Prefs.get("com.braincadet.phd.ni", ni_csv), 10);
            gd.addStringField("step", Prefs.get("com.braincadet.phd.step", step_csv), 10);
            gd.addStringField("kappa", Prefs.get("com.braincadet.phd.kappa", kappa_csv), 10);
            gd.addStringField("ps", Prefs.get("com.braincadet.phd.ps", pS_csv), 10);
            gd.addStringField("pd", Prefs.get("com.braincadet.phd.pd", pD_csv), 10);
            gd.addStringField("krad", Prefs.get("com.braincadet.phd.krad", krad_csv), 10);
            gd.addStringField("kc", Prefs.get("com.braincadet.phd.kc", kc_csv), 10);
            gd.addNumericField("maxiter", Prefs.get("com.braincadet.phd.maxiter", maxiter), 0, 5, "");
            gd.addStringField("maxepoch", Prefs.get("com.braincadet.phd.maxepoch", Integer.toString(maxepoch)), 10);
            gd.addCheckbox("savemidres", Prefs.get("com.braincadet.phd.savemidres", savemidres));
            gd.addCheckbox("usetness", Prefs.get("com.braincadet.phd.usetness", usetness));

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
            maxepoch = (maxepoch_str.equals("Inf")) ? Integer.MAX_VALUE : Integer.valueOf(maxepoch_str);
            Prefs.set("com.braincadet.phd.maxepoch", maxepoch);
            savemidres = gd.getNextBoolean();
            Prefs.set("com.braincadet.phd.savemidres", savemidres);
            usetness = gd.getNextBoolean();
            Prefs.set("com.braincadet.phd.usetness", usetness);

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
            maxepoch = (maxepoch_str.equals("Inf")) ? Integer.MAX_VALUE : Integer.valueOf(maxepoch_str);
            savemidres = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "savemidres", String.valueOf(false)));
            usetness = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "usetness", String.valueOf(true)));
        }

        if (savemidres) {
            midresdir = ip_load.getOriginalFileInfo().directory + ip_load.getTitle() + "_midres";//File.separator +
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

//        dd = maxepoch_csv.split(","); if (dd.length==0) return;
//        maxepoch = new int[dd.length];
//        for (int i = 0; i < dd.length; i++) maxepoch[i] = Integer.valueOf(dd[i]);
//        if (true) {IJ.log("ok, maxepoch = " + maxepoch); return;}

        //******************************************************************
        ImageStack is_tness;
        ImagePlus ip_tness;

        IJ.log(" -- prefiltering...");

        long t1prep = System.currentTimeMillis();

        is_tness = new ImageStack(N, M);
        for (int i = 0; i < P; i++) {
            float[] tt = new float[N * M];
            Arrays.fill(tt, 0f);
            is_tness.addSlice(new FloatProcessor(N, M, tt));
        }

        ip_tness = new ImagePlus("tness", is_tness);
        String[] readLn = sigmas.trim().split(",");

        for (int i = 0; i < readLn.length; i++) {
            float sig = Float.valueOf(readLn[i].trim()).floatValue();
            IJ.log("sig=" + IJ.d2s(sig,2));
            TubenessProcessor tp = new TubenessProcessor(sig, false);
            ImagePlus result = tp.generateImage(ip_load);
            ImageCalculator ic = new ImageCalculator();
            // average
//          IJ.run(result, "Multiply...", "value=" + IJ.d2s(1f/readLn.length,3) + " stack");
//          ic.run("Add 32-bit stack", ip_tness, result); // result of the addition is placed in ip_tness
            // max
            ic.run("Max 32-bit stack", ip_tness, result);
        }

        ip_tness.setCalibration(null);

        if (savemidres) {
            ImagePlus temp = ip_tness.duplicate();
            IJ.run(temp, "8-bit", ""); // convert to 8 bit before saving
//                IJ.log("saving... " + midresdir + File.separator + "tness," + sigmas + ".tif");
            IJ.saveAs(temp, "Tiff", midresdir + File.separator + "tness," + sigmas + ".tif");
        }

        // tubeness min-max normalize and store in an array for later and extract locations in a separate array
        float tnessmin = Float.POSITIVE_INFINITY;
        float tnessmax = Float.NEGATIVE_INFINITY;
        tness = new float[SZ];
        int[][] locationXYZ = new int[SZ][3]; // random sampling weighted with the normalized tness as importance function

        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

            float[] slc_float = null;
            byte[] slc_byte = null;

            if (usetness)
                slc_float = (float[]) ip_tness.getStack().getPixels(z);
            else
                slc_byte = (byte[]) ip_load.getStack().getPixels(z);

            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {
                    int ii = (z - 1) * (N * M) + y * N + x;

                    if (usetness)
                        tness[ii] = slc_float[y * N + x];
                    else
                        tness[ii] = (float) (slc_byte[y * N + x] & 0xff);

                    locationXYZ[ii][0] = x;
                    locationXYZ[ii][1] = y;
                    locationXYZ[ii][2] = (z - 1);

                    if (tness[ii] < tnessmin) tnessmin = tness[ii];
                    if (tness[ii] > tnessmax) tnessmax = tness[ii];

                }
            }
        }

        for (int i = 0; i < SZ; i++) {
            tness[i] = (tnessmax - tnessmin > Float.MIN_VALUE) ? ((tness[i] - tnessmin) / (tnessmax - tnessmin)) : 0;
        }

        long t2prep = System.currentTimeMillis();
        IJ.log("t_preprocessing" + IJ.d2s((t2prep - t1prep) / 1000f,2) + "[sec]");

        suppmap = new int[SZ];         // suppression map with node tags

        // go through comma separated parameter values
        for (int i01 = 0; i01 < no.length; i01++) {
            for (int i02 = 0; i02 < ro.length; i02++) {
                for (int i03 = 0; i03 < ni.length; i03++) {
                    for (int i04 = 0; i04 < krad.length; i04++) {
                        for (int i05 = 0; i05 < step.length; i05++) {
                            for (int i06 = 0; i06 < kappa.length; i06++) {
                                for (int i07 = 0; i07 < pS.length; i07++) {
                                    for (int i08 = 0; i08 < pD.length; i08++) {
                                        for (int i09 = 0; i09 < th.length; i09++) { // todo: change threshold into local maxima sensitivity
                                            for (int i10 = 0; i10 < kc.length; i10++) {
//                                                for (int i11 = 0; i11 < maxepoch.length; i11++) {

                                                long t1 = System.currentTimeMillis();

                                                if (savemidres) {
                                                    Tools.createAndCleanDir(midresdir + File.separator + "g(z|x)");
                                                    Tools.createAndCleanDir(midresdir + File.separator + "suppmap");
                                                    Tools.createAndCleanDir(midresdir + File.separator + "objects");

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

//                                                    ImagePlus impool = ip_tness.duplicate();//(usetness)?ip_tness.duplicate():ip_load.duplicate();
//                                                    IJ.run(impool, "8-bit", "");

//                                                    int threshold = (int) Math.ceil(th[i09] * 255);
//                                                    applythreshold(threshold, impool);

//                                                    Prefs.blackBackground = true;
//                                                    IJ.run(impool, "Skeletonize", "stack");

//                                                    if (savemidres) {
//                                                        IJ.saveAs(impool, "Tiff", midresdir + File.separator + "seedpool,th=" + IJ.d2s(th[i09], 2) + ".tif");
//                                                    }
                                                ArrayList<Integer> locs = new ArrayList<Integer>();     // list of candidate locations for seed points
                                                ArrayList<Float> locsw = new ArrayList<Float>();       // weights assigned to each location

                                                for (int z = ((P == 1) ? 1 : 2); z <= ((P == 1) ? P : P - 1); z++) { // layer count, zcoord is layer-1
//                                                        byte[] slc;
                                                    Polygon maxx;
//                                                        if (true) { // seed location - local maxima of the tubularity measure
                                                    MaximumFinder mf = new MaximumFinder();
//                                                            slc = (byte[])mf.findMaxima(ip_tness.getStack().getProcessor(z), 0.1, MaximumFinder.SINGLE_POINTS, true).getPixels();
                                                    maxx = mf.getMaxima(ip_tness.getStack().getProcessor(z), th[i09], false);
//                                                        }
//                                                        else { // alternative seed location - threshold+skeletonize tubularity measure
//                                                            slc = (byte[]) impool.getStack().getPixels(z);
//                                                        }

                                                    for (int i = 0; i < maxx.npoints; i++) {
//                                                            IJ.log("[x,y]="+maxx.xpoints[i]+", "+maxx.ypoints[i]);
                                                        int ii = (z - 1) * (N * M) + maxx.ypoints[i] * N + maxx.xpoints[i];
                                                        locs.add(ii);
                                                        locsw.add((float) Math.pow(tness[ii], MultiTT.weight_deg)); // 1f (float) Math.pow(tness[ii], MultiTT.weight_deg77)
                                                    }
//                                                        for (int x = 0; x < N; x++) {
//                                                            for (int y = 0; y < M; y++) {
//                                                                int ii = (z - 1) * (N * M) + y * N + x;
//                                                                if ((slc[y * N + x] & 0xff) == 255) {
//                                                                    locs.add(ii);
//                                                                    locsw.add((float) Math.pow(tness[ii], MultiTT.weight_deg77));
//                                                                }
//                                                            }
//                                                        }
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

//                                                    Overlay ov = new Overlay();
//                                                    OvalRoi or = new OvalRoi(5 + 0.5, 5 + 0.5, 2, 2);
//                                                    or.setFillColor(Color.YELLOW);
//                                                    or.setPosition(0 + 1);
//                                                    ov.add(or);
                                                    // save input image with an overlay of candidate points
//                                                        ip_load.setOverlay(ov);
//                                                        IJ.saveAs(ip_load, "TIFF", midresdir + File.separator + "t0.tif");

                                                    exportlocsxyz(tt, 0.3f, VIOLET, midresdir, "seedpool_tolerance=" + IJ.d2s(th[i09], 2));

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

                                                // multi-object detection video demo
//                                                template_stack = new ImageStack(N, M);
//                                                for (int lay = 0; lay < P; lay++) {
//                                                    template_stack.addSlice(new ByteProcessor(N, M));
//                                                }
//                                                template_image = new ImagePlus();
//                                                template_ovrly = new Overlay();
//                                                template_zprojector = new ZProjector();
                                                //******************************************************************
                                                IJ.log("-- multi-object filtering...");
//                                                long t1 = System.currentTimeMillis();

                                                int epochcnt = 0;

                                                while (locs.size() > 0 && epochcnt < maxepoch) {

                                                    epochcnt++;

                                                    if (mtt.verbose) IJ.log("e=" + epochcnt + " [" + maxepoch + "]");

                                                    iter_count = 0;

                                                    int cnt_removed = 0;

                                                    for (int i = locs.size() - 1; i >= 0; i--) { // exclude filtering from covered voxels
                                                        if (suppmap[locs.get(i)] > 0) {
                                                            locs.remove(i);
                                                            locsw.remove(i);
                                                            cnt_removed++;
                                                        }
                                                    }

                                                    if (mtt.verbose)
                                                        IJ.log(IJ.d2s(locs.size() / 1000f, 1) + "k locations [" + locs.size() + "] \n" + IJ.d2s((locs.size() / locs_count) * 100f, 1) + "% of the initial pool\n------------------\n");

                                                    if (locs.size() == 0) {
                                                        IJ.log("locs.size()==0");

                                                        // export tree
                                                        String delindir = imdir + "PHD.sig.th.no.ro.ni.krad.stp.kapa.ps.pd.kc.e_" + sigmas + "_" + IJ.d2s(th[i09], 2) + "_" + IJ.d2s(no[i01], 0) + "_" + IJ.d2s(ro[i02], 0) + "_" + IJ.d2s(ni[i03], 0) + "_" + IJ.d2s(krad[i04], 0) + "_" + IJ.d2s(step[i05], 0) + "_" + IJ.d2s(kappa[i06], 1) + "_" + IJ.d2s(pS[i07], 2) + "_" + IJ.d2s(pD[i08], 2) + "_" + IJ.d2s(kc[i10], 1) + "_" + IJ.d2s(epochcnt, 0);// + "_" + IJ.d2s(new Random().nextInt(Integer.MAX_VALUE),0);
//                                                            Tools.createAndCleanDir(delindir);
                                                        Tools.createDir(delindir);
                                                        remove_double_links(mtt.Y);
                                                        if (!is_biderectinal_linking(mtt.Y)) {
                                                            IJ.log("missing link! fault in bidirectional linking");
                                                            return;
                                                        }

                                                        ArrayList<Node> tree = mtt.bfs1(mtt.Y, true);
                                                        exportReconstruction(tree, delindir, imnameshort); // export .swc   epochcnt, mtt.Y.size()-1,

                                                        break; // go out of while()
                                                    }

//                                                        ArrayList<int[]> N_o = initlocs(no[i01], locs, locsw); // xyz locations
                                                    ArrayList<int[]> N_o = new ArrayList<int[]>();

                                                    //********** initialization **********//
                                                    mtt._init(no[i01], step[i05], locs, locsw, N_o, img, N, M, P, tness, suppmap); // , template_ovrly
                                                    if (N_o.size() == 0) {
                                                        IJ.log("initialization stopped, |N_o|=0");
                                                        break; // out of while()
                                                    }

                                                    //exportlocsxyz(N_o, 10f, RED, ip_load.getOriginalFileInfo().directory + File.separator, IJ.d2s(N_o.size(),0)+"_seeds_");

                                                    if (savemidres) {

                                                        exportlocsxyz(N_o, 5f, RED, midresdir, "seeds,e=" + IJ.d2s(epochcnt, 0));

                                                        exportXYZW(mtt.Xk, midresdir, "XWinit,epoch=" + IJ.d2s(epochcnt, 0), MAGENTA);      // phd weights
                                                        exportXYZVxyz(mtt.Xk, midresdir, "XVinit,epoch=" + IJ.d2s(epochcnt, 0), BLUE);     // directions of particles
                                                        exportNodes(mtt.Y, midresdir, "Yinit,epoch=" + IJ.d2s(epochcnt, 0)); // estimations

                                                        logval(tnessCsvLog, mtt.cluttertness);
                                                        logval(tnessCsvLog, mtt.Xk);
                                                        logval(zsizeCsvLog, no[i01]);                     // nr. observations mtt.Zk.size()
                                                        logval(phdmassCsvLog, mtt.phdmass);

                                                        //---------------------------------------------------------------------------------------------------
                                                        // multi-object detection video demo
//                                                        if (false) { // it can log the tags... slows down a lot... and takes memory!!!
//                                                            ImagePlus hpimp = getSuppMap();
//                                                            template_zprojector.setImage(hpimp);
//                                                            template_zprojector.setMethod(ZProjector.MAX_METHOD);
//                                                            template_zprojector.doProjection();
//                                                            hpimp = template_zprojector.getProjection();
//                                                                IJ.run(hpimp, "8-bit", "");
//                                                            IJ.saveAs(hpimp, "Zip", midresdir + File.separator + "suppmap" + File.separator + "r=" + IJ.d2s(epochcnt, 0) + ",i=" + IJ.d2s(0, 0) + ".zip");
//                                                                template_image.setStack(template_stack);
//                                                                template_image.setOverlay(template_ovrly);
//                                                                template_image.flattenStack(); // takes java 1.6 at least
//                                                                IJ.run(template_image, "8-bit", "");
//                                                                template_zprojector.setImage(template_image);
//                                                                template_zprojector.setMethod(ZProjector.MAX_METHOD);
//                                                                template_zprojector.doRGBProjection();
//                                                                ImagePlus hpimp1 = template_zprojector.getProjection();

//                                                                hpimp1.setTitle("obj");
//                                                                hpimp1.show();
//                                                                hpimp.show();
//                                                                IJ.run("Add Image...", "image=obj x=0 y=0 opacity=40 zero");
//                                                                hpimp.updateAndDraw();

//                                                                IJ.saveAs(hpimp, "Zip", midresdir + File.separator + "suppmap" + File.separator + "phd,r="+IJ.d2s(epochcnt, 0)+",i=" + IJ.d2s(0, 0) + ".zip");

//                                                                hpimp.close();
//                                                                hpimp1.close();
//                                                        }

                                                    }

                                                    if (mtt.Xk.size() > 0) {

//                                                            IJ.log("mtt.Xk.size() > 0");
//                                                            IJ.log("mtt.Y.size()="+mtt.Y.size());
//                                                            if (savemidres) Xlog(mtt.Xk, GREEN);

                                                        while (iter_count < maxiter) {

                                                            if (mtt.verbose) IJ.log("k=" + iter_count);

                                                            boolean iterok;

                                                            iterok = mtt._iter1(N, M, P, tness, suppmap); // multi-object detection video demo: , template_ovrly

                                                            if (savemidres) { // iter_count%5==0 && iter_count==maxiter-1

                                                                Xlog(mtt.Xk, GREEN);
                                                                XPlog(mtt.XPk, BLUEBERRY);
                                                                Zlog(mtt.Zk, RED, 1f);
                                                                ZPlog(mtt.ZPk, OCHRE, .1f);

                                                                //---------------------------------------------------------------------------------------------------
                                                                // multi-object detection video demo
//                                                                if (false) { // set if you wish to have the map, takes lot of resources!

//                                                                    ImagePlus hpimp = getSuppMap();
//                                                                    template_zprojector.setImage(hpimp);
//                                                                    template_zprojector.setMethod(ZProjector.MAX_METHOD);
//                                                                    template_zprojector.doProjection();
//                                                                    hpimp = template_zprojector.getProjection();
//                                                                        IJ.run(hpimp, "8-bit", "");
//                                                                    IJ.saveAs(hpimp, "Zip", midresdir + File.separator + "suppmap" + File.separator + "r=" + IJ.d2s(epochcnt, 0) + ",i=" + IJ.d2s(iter_count + 1, 0) + ".zip");

//                                                                        template_image.setStack(template_stack);
//                                                                        if (template_ovrly.size()>0) {
//                                                                            template_image.setOverlay(template_ovrly);
//                                                                            template_image.flattenStack(); // asks java 1.6
//                                                                        }

//                                                                      IJ.run(template_image, "8-bit", "");
//                                                                        template_zprojector.setImage(template_image);
//                                                                        template_zprojector.setMethod(ZProjector.MAX_METHOD);
//                                                                        if (template_image.getType()==ImagePlus.COLOR_RGB)
//                                                                            template_zprojector.doRGBProjection(); // there was flattening
//                                                                        else
//                                                                            template_zprojector.doProjection();

//                                                                        ImagePlus hpimp1 = template_zprojector.getProjection();

//                                                                        hpimp1.setTitle("obj");
//                                                                        hpimp1.show();
//                                                                        hpimp.show();
//                                                                        IJ.run("Add Image...", "image=obj x=0 y=0 opacity=40 zero");
//                                                                        hpimp.updateAndDraw();

//                                                                        IJ.saveAs(hpimp, "Zip",midresdir + File.separator + "suppmap" + File.separator + "phd,r="+IJ.d2s(epochcnt, 0)+",i=" + IJ.d2s(iter_count+1, 0) + ".zip");

//                                                                        hpimp.close();
//                                                                        hpimp1.close();

//                                                                }

                                                                logval(tnessCsvLog, mtt.XPk);
                                                                logval(zsizeCsvLog, mtt.Zk.size());
                                                                logval(phdmassCsvLog, mtt.phdmass);

                                                                if (mtt.g != null) {
                                                                    ImagePlus gimp = new ImagePlus("g(z|x),iter0=" + IJ.d2s(iter_count, 0), new FloatProcessor(mtt.g));
                                                                    IJ.run(gimp, "Rotate 90 Degrees Right", "");
                                                                    IJ.saveAs(gimp, "Tiff", midresdir + File.separator + "g(z|x)" + File.separator + "g(z|x),iter0=" + IJ.d2s(iter_count, 0) + ".tif");
                                                                }

                                                            }

                                                            if (!iterok)
                                                                break; // go out of while loop and try new set of seed points

                                                            iter_count++;

                                                        }
                                                    } else IJ.log("mtt.Xk.size() == 0");

                                                    if (locs.size() == 0 || epochcnt == maxepoch) { // each # of iterations  || epochcnt%5 == 0

                                                        IJ.log("export the nodes....");

                                                        // dev.
                                                        IJ.log("before: mtt.Y.size =" + mtt.Y.size());

                                                        interpolate_nodelist(mtt.Y, 1f);

                                                        IJ.log("after: mtt.Y.size = " + mtt.Y.size());



                                                        if (true) return;

                                                        // save output
                                                        String delindir = imdir + "PHD.sig.th.no.ro.ni.krad.stp.kapa.ps.pd.kc.e_" + sigmas + "_" + IJ.d2s(th[i09], 2) + "_" + IJ.d2s(no[i01], 0) + "_" + IJ.d2s(ro[i02], 0) + "_" + IJ.d2s(ni[i03], 0) + "_" + IJ.d2s(krad[i04], 0) + "_" + IJ.d2s(step[i05], 0) + "_" + IJ.d2s(kappa[i06], 1) + "_" + IJ.d2s(pS[i07], 2) + "_" + IJ.d2s(pD[i08], 2) + "_" + IJ.d2s(kc[i10], 1) + "_" + IJ.d2s(epochcnt, 0);// + "_" + IJ.d2s(new Random().nextInt(Integer.MAX_VALUE),0);
                                                        Tools.createDir(delindir);
                                                        remove_double_links(mtt.Y);
                                                        if (!is_biderectinal_linking(mtt.Y)) {
                                                            IJ.log("missing link! fault in bidirectional linking");
                                                            return;
                                                        }

                                                        ArrayList<Node> tree = mtt.bfs1(mtt.Y, true);
                                                        exportReconstruction(tree, delindir, imnameshort); // export .swc   epochcnt, mtt.Y.size()-1,

                                                        if (savemidres) {
                                                            ImagePlus hpimp = getSuppMap();
                                                            IJ.saveAs(hpimp, "Zip", midresdir + File.separator + "suppmap,rounds=" + IJ.d2s(epochcnt, 0) + ",i=" + IJ.d2s(iter_count + 1, 0) + ".zip"); // "suppmap" + File.separator +
                                                        }
                                                    }

                                                } // while there are locations and epochs have not reached the limit

                                                long t2 = System.currentTimeMillis();
                                                IJ.log("done. " + IJ.d2s(((t2 - t1) / 1000f), 2) + "s. [maxiter=" + maxiter + ", rounds=" + maxepoch + "]");

                                                if (savemidres) {
                                                    exportDelineation(mtt.Y, midresdir, imnameshort);   // .phd file
                                                    exportNodes(mtt.Y, midresdir, imnameshort);         // .swc file with isolated nodes
                                                }

                                                // clear mtt components
                                                mtt.Xk.clear();
                                                mtt.Y.clear();

//                                                }
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

    private void interpolate_nodelist(ArrayList<Node> nlist, float step) {

        // assume bidirectional connections in nlist

        // interpolate all inter-node links with the step size
        ArrayList< ArrayList<Boolean> > chk = new ArrayList<>(2*nlist.size());
        for (int i = 0; i < nlist.size(); i++) {
            ArrayList<Boolean> chk1 = new ArrayList<>();
            if (nlist.get(i)!=null) {
                for (int j = 0; j < nlist.get(i).nbr.size(); j++) {
                    chk1.add(false);
                }
            }
            chk.add(chk1);
        }

        int init_size = nlist.size();

        for (int i = 1; i < init_size; i++) {
            for (int j = 0; j < nlist.get(i).nbr.size(); j++) {
                if (!chk.get(i).get(j)) {

                    int i1 = nlist.get(i).nbr.get(j);
                    int j1 = nlist.get(i1).nbr.indexOf(i); // find(nX[i1].nbr.begin(), nX[i1].nbr.end(), i) - nX[i1].nbr.begin();

                    if (j1 != -1) { // interpolate if there was existing link back

                        chk.get(i).set(j, true);
                        chk.get(i1).set(j1, true);

                        float vnorm = (float) Math.sqrt(Math.pow(nlist.get(i1).loc[0]-nlist.get(i).loc[0], 2) +
                                Math.pow(nlist.get(i1).loc[1]-nlist.get(i).loc[1], 2) +
                                Math.pow(nlist.get(i1).loc[2]-nlist.get(i).loc[2], 2)
                        );

                        float vx = (nlist.get(i1).loc[0]-nlist.get(i).loc[0])/vnorm;
                        float vy = (nlist.get(i1).loc[1]-nlist.get(i).loc[1])/vnorm;
                        float vz = (nlist.get(i1).loc[2]-nlist.get(i).loc[2])/vnorm;
                        int N = (int) Math.ceil(vnorm/step);

                        // add subsampling if N>1
                        for (int k = 1; k < N; ++k) {
                            // add the node,only location is used in the refinement stage currently
                            nlist.add(new Node(
                                    nlist.get(i).loc[0] + k * (vnorm/N) * vx,
                                    nlist.get(i).loc[1] + k * (vnorm/N) * vy,
                                    nlist.get(i).loc[2] + k * (vnorm/N) * vz,
                                    nlist.get(i).r + (nlist.get(i1).r - nlist.get(i).r) * (k/(float)N),
                                    ((k<=N/2)?nlist.get(i).type:nlist.get(i1).type)
                            ));

                            // link backward
                            if (k==1) {
                                // first: link nlist[nlist.size()-1] with nlist[i]
                                nlist.get(nlist.size()-1).nbr.add(i); // nX[nX.size()-1].nbr.push_back(i);
                                nlist.get(i).nbr.set(j, nlist.size()-1); //  nX[i].nbr[j] = nX.size()-1; // replace i1 with the link to the first addition
                            }
                            else {
                                // middle: link nlist[nlist.size()-1] with nlist[nlist.size()-2]
                                nlist.get(nlist.size()-1).nbr.add(nlist.size()-2); // nX[nX.size()-1].nbr.push_back(nX.size()-2);
                                nlist.get(nlist.size()-2).nbr.add(nlist.size()-1); // nX[nX.size()-2].nbr.push_back(nX.size()-1);
                            }

                            // link forward
                            if (k==N-1) {
                                // last: link nlist[nlist.size()-1] with nlist[i1]
                                nlist.get(nlist.size()-1).nbr.add(i1); // nX[nX.size()-1].nbr.push_back(i1);
                                nlist.get(i1).nbr.set(j1, nlist.size()-1); // nX[i1].nbr[j1] = nX.size()-1; // replace i with the link to the last addition
                            }

                        }

                    }

                }
            }
        }

        IJ.log("" + IJ.d2s(((float)nlist.size()/init_size)*100f, 2) + " % node # after interpolation");

    }

    private void interpolate_treelist(ArrayList<Node> ntree, float step, int type) {

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
            for (int j = 0; j < ntree[i].nbr.size(); ++j) { // there should be 0 or 1 neighbor
                int i1 = ntree.get(i).nbr.get(j);

                // interpolate between ntree[i] and ntree[i1]
                float vnorm = (float) Math.sqrt(Math.pow(ntree.get(i1).loc[0]-ntree.get(i).loc[0], 2) +
                        Math.pow(ntree.get(i1).loc[1]-ntree.get(i).loc[1], 2) +
                        Math.pow(ntree.get(i1).loc[2]-ntree.get(i).loc[2], 2));
                float vx = (ntree.get(i1).loc[0]-ntree.get(i).loc[0])/vnorm;
                float vy = (ntree.get(i1).loc[1]-ntree.get(i).loc[1])/vnorm;
                float vz = (ntree.get(i1).loc[2]-ntree.get(i).loc[2])/vnorm;
                int N = (int) Math.ceil(vnorm/step);

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

    private void non_blurring(ArrayList<Node> nX, ArrayList<Node> nY, float SIG2RAD, int MAXITER, float EPSILON2) {

        // mean-shift (non-blurring) uses flexible neighbourhood scaled with respect to the node's sigma
        int checkpoint = (int) Math.round(nX.size()/10.0);

        float[] conv = {0,0,0,0}; // x y z r
        float[] next = {0,0,0,0}; // x y z r


    }

    private void remove_double_links(ArrayList<Node> nlist) {
        for (int i = 0; i < nlist.size(); i++) {
            if (nlist.get(i) != null) {
                Set<Integer> set = new HashSet<Integer>();
                set.addAll(nlist.get(i).nbr);
                nlist.get(i).nbr.clear();
                nlist.get(i).nbr.addAll(set);
            }
        }
    }

    private boolean is_biderectinal_linking(ArrayList<Node> nlist) {
        for (int i = 0; i < nlist.size(); i++) {
            if (nlist.get(i) != null) {
                for (int j = 0; j < nlist.get(i).nbr.size(); j++) {
                    int nbr_idx = nlist.get(i).nbr.get(j);
                    if (Collections.frequency(nlist.get(nbr_idx).nbr, i) != 1) {
                        IJ.log("ERROR: " + i + " --> " + nbr_idx);
                        return false;
                    }
                }
            }
        }
        return true;
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

        } catch (IOException e) {
        }

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

    // for the publication report only
//    private void exportTime(String name, float t22, String signature, String outdir, String outfile) {}

    private void exportDelineation(ArrayList<Node> nlist, String outdir, String outfile) {

        Tools.createDir(outdir);
        String delinswc1 = outdir + File.separator + outfile + ".phd";
        Tools.cleanfile(delinswc1);

        try {

            PrintWriter logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(delinswc1, true)));

            int t1 = 0;

            for (int i = 0; i < nlist.size(); i++) {
                if (nlist.get(i) != null) {

                    logWriter1.print((++t1) + " " + YELLOW + " " + IJ.d2s(nlist.get(i).loc[0], 3) + " " + IJ.d2s(nlist.get(i).loc[1], 3) + " " + IJ.d2s(nlist.get(i).loc[2], 3) + " " + IJ.d2s(nlist.get(i).r, 3) + " ");
                    for (int j = 0; j < nlist.get(i).nbr.size(); j++) {
                        logWriter1.print(IJ.d2s(nlist.get(i).nbr.get(j), 0) + " ");
                    }
                    logWriter1.println("");

                }
            }

            logWriter1.close();

        } catch (IOException e) {
        }

    }

    private void exportReconstruction(ArrayList<Node> nlist, String outdir, String outname) {

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

    public ImagePlus getSuppMap() { // String title

        ImageStack outis = new ImageStack(N, M);

        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

            float[] slc = new float[N * M];

            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {

                    int ii = (z - 1) * (N * M) + y * N + x;

                    slc[y * N + x] = suppmap[ii];

                }
            }

            outis.addSlice(new FloatProcessor(N, M, slc));

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

}
