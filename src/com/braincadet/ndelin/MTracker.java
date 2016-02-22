package com.braincadet.ndelin;

import com.braincadet.ndelin.fun.Tools;
import com.braincadet.ndelin.multi.MultiTT;
import com.braincadet.ndelin.multi.X;
import com.braincadet.ndelin.swc.Node;
import features.TubenessProcessor;
import ij.*;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.measure.Measurements;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.*;

import java.io.*;
import java.util.*;

public class MTracker implements PlugIn {

    // color coding for the swc in vaa3d
    static int BLUE=3, YELLOW=6, RED=2, BLACK=1, MAGENTA=5, WHITE=0, VIOLET=4, GREEN=7, OCHRE=8, WEAK_GREEN=9, PINK=10, BLUEBERRY=12;

    AutoThresholder.Method thmethod = AutoThresholder.Method.IsoData; // IsoData Moments Otsu Triangle Default

    // allow having array of comma separated parameter inputs for the tubularity measure
    String sigmas = "2,3";                   // comma separated scale values (tubularity)

    // one sigmas tubularity image can be called for sequence of parameters (comma separated values)
    // parameter lists (string + array with extracted values)
    int[] no;                           // initial multi-object state cardinality
    String no_csv = "";                 // comma separated string of test parameters

    int[]   ro;                  // number of particles per object
    String  ro_csv = "";

    int[]   ni;                  // number of predictions per particle
    String  ni_csv = "";

    int[] krad;                  // tube diameter in pixels
    String krad_csv = "";

    int[]   step;                // motion step
    String  step_csv = "";

    float[] kappa;                      // von Mises angular probability
    String  kappa_csv = "";

    float[] pS;                         // survival pty
    String  pS_csv = "";

    float[] pD;                         // detection pty
    String  pD_csv = "";

    float[] th;                         // reference tubularity ratio for clutter
    String th_csv = "";

    float[] kc;                         // clutter phd decay parameter
    String kc_csv = "";

    int[]  maxepoch;                    // epoch limit
    String maxepoch_csv = "";           // comma separated input

    int TUBE_RADIUS = 10;                // used at _init(), radius of the sphere used for the seed point

    // save results
    int     maxiter = -1;               // iteration limit (hundreds are fine)
    boolean savemidres = false;         // save partial results
    boolean usetness = true;            // use tubularity measure

    float[] img;                        // original image as float array
    float[] tness;                      // tubeness min-max normalized
    int[]   suppmap;                    // supression map: disable sampling (image stack size)

    int N, M, P, SZ;                    // stack dimensions (width, height, length, size)
    String imdir, imnameshort;
    String midresdir = "";              // output directories, filenames

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

    private static final long MEGABYTE = 1024L * 1024L;

    public static long bytesToMegabytes(long bytes) {
        return bytes / MEGABYTE;
    }

    public void run(String s) {

        // read input image
        String in_folder = Prefs.get("com.braincadet.ndelin.multi.dir", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory(); Prefs.set("com.braincadet.ndelin.multi.dir", in_folder);
        String image_path = dc.getPath();
        if (image_path==null) return;

        ImagePlus ip_load = new ImagePlus(image_path);

        if(ip_load==null) {
            IJ.log(image_path + " was null");
            return;
        }
        if (ip_load.getType()!=ImagePlus.GRAY8) {
            IJ.log("Image needs to be GRAY8.");
            return;
        }

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

        if (Macro.getOptions()==null) {

            GenericDialog gd = new GenericDialog("TreeDelinPHD");
            gd.addMessage("");//tubularity measure
            gd.addStringField("sigmas",       Prefs.get("com.braincadet.ndelin.multi.sigmas",       sigmas), 10);
            gd.addStringField("th",           Prefs.get("com.braincadet.ndelin.multi.th",           th_csv), 10);
            gd.addMessage("");
            gd.addStringField("no",           Prefs.get("com.braincadet.ndelin.multi.no",           no_csv), 10);
            gd.addStringField("ro",           Prefs.get("com.braincadet.ndelin.multi.ro",           ro_csv), 10);
            gd.addStringField("ni",             Prefs.get("com.braincadet.ndelin.multi.ni",           ni_csv), 10);
            gd.addStringField("step",           Prefs.get("com.braincadet.ndelin.multi.step",         step_csv), 10);
            gd.addStringField("kappa",          Prefs.get("com.braincadet.ndelin.multi.kappa",        kappa_csv), 10);
            gd.addStringField("ps",             Prefs.get("com.braincadet.ndelin.multi.ps",           pS_csv), 10);
            gd.addStringField("pd",             Prefs.get("com.braincadet.ndelin.multi.pd",           pD_csv), 10);
            gd.addMessage("");//clustering kernel
            gd.addStringField("krad",           Prefs.get("com.braincadet.ndelin.multi.krad",         krad_csv), 10);
            gd.addMessage("");//clutter
            gd.addStringField("kc",             Prefs.get("com.braincadet.ndelin.multi.kc",           kc_csv), 10);
            gd.addMessage("");
            gd.addNumericField("maxiter",       Prefs.get("com.braincadet.ndelin.multi.maxiter",      maxiter), 0, 5, "");
            gd.addStringField("maxepoch",       Prefs.get("com.braincadet.ndelin.multi.maxepoch",     maxepoch_csv), 10);
            gd.addCheckbox("savemidres",        Prefs.get("com.braincadet.ndelin.multi.savemidres",   savemidres));
            gd.addCheckbox("usetness",          Prefs.get("com.braincadet.ndelin.multi.usetness",     usetness));

            gd.showDialog();
            if (gd.wasCanceled()) return;

            sigmas = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.sigmas",     sigmas);
            th_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.th",         th_csv);
            no_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.no",         no_csv);
            ro_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.ro",         ro_csv);
            ni_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.ni",         ni_csv);
            step_csv = gd.getNextString();        Prefs.set("com.braincadet.ndelin.multi.step",       step_csv);
            kappa_csv = gd.getNextString();       Prefs.set("com.braincadet.ndelin.multi.kappa",      kappa_csv);
            pS_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.ps",         pS_csv);
            pD_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.pd",         pD_csv);
            krad_csv = gd.getNextString();        Prefs.set("com.braincadet.ndelin.multi.krad",       krad_csv);
            kc_csv = gd.getNextString();          Prefs.set("com.braincadet.ndelin.multi.kc",         kc_csv);
            maxiter = (int) gd.getNextNumber();   Prefs.set("com.braincadet.ndelin.multi.maxiter",    maxiter);
            maxepoch_csv = gd.getNextString();    Prefs.set("com.braincadet.ndelin.multi.maxepoch",   maxepoch_csv);
            savemidres = gd.getNextBoolean();     Prefs.set("com.braincadet.ndelin.multi.savemidres", savemidres);
            usetness = gd.getNextBoolean();       Prefs.set("com.braincadet.ndelin.multi.usetness",   usetness);

        }
        else {

            sigmas = Macro.getValue(Macro.getOptions(), "sigmas",                                     sigmas);
            th_csv = Macro.getValue(Macro.getOptions(), "th",                                         th_csv);
            no_csv = Macro.getValue(Macro.getOptions(), "no",                                         no_csv);
            ro_csv = Macro.getValue(Macro.getOptions(), "ro",                                         ro_csv);
            ni_csv = Macro.getValue(Macro.getOptions(), "ni",                                         ni_csv);
            step_csv = Macro.getValue(Macro.getOptions(), "step",                                     step_csv);
            kappa_csv = Macro.getValue(Macro.getOptions(), "kappa",                                   kappa_csv);
            pS_csv = Macro.getValue(Macro.getOptions(), "ps",                                         pS_csv);
            pD_csv = Macro.getValue(Macro.getOptions(), "pd",                                         pD_csv);

            krad_csv = Macro.getValue(Macro.getOptions(), "krad",                                     krad_csv);

            kc_csv = Macro.getValue(Macro.getOptions(), "kc",                                         kc_csv);

            maxiter = Integer.valueOf(Macro.getValue(Macro.getOptions(),     "maxiter",               String.valueOf(maxiter)));
            maxepoch_csv = Macro.getValue(Macro.getOptions(), "maxepoch",                             maxepoch_csv);
            savemidres = Boolean.valueOf(Macro.getValue(Macro.getOptions(),  "savemidres",            String.valueOf(false)));
            usetness = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "usetness",                 String.valueOf(true)));
        }

        if (savemidres) {
            midresdir = ip_load.getOriginalFileInfo().directory + ip_load.getTitle() + "_midres";//File.separator +
            Tools.createAndCleanDir(midresdir); // create midresult dir and initialize export/log
        }

        String[] dd;

        // ge comma separated parameter values
        dd = th_csv.split(","); if (dd.length==0) return;
        th = new float[dd.length];
        for (int i = 0; i < dd.length; i++) th[i] = Float.valueOf(dd[i]);

        dd = no_csv.split(",");  if (dd.length==0) return;
        no = new int[dd.length];
        for (int i = 0; i < dd.length; i++) no[i] = Integer.valueOf(dd[i]);

        dd = ro_csv.split(","); if (dd.length==0) return;
        ro = new int[dd.length];
        for (int i = 0; i < dd.length; i++) ro[i] = Integer.valueOf(dd[i]);

        dd = ni_csv.split(","); if (dd.length==0) return;
        ni = new int[dd.length];
        for (int i = 0; i < dd.length; i++) ni[i] = Integer.valueOf(dd[i]);

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

        dd = krad_csv.split(","); if (dd.length==0) return;
        krad = new int[dd.length];
        for (int i = 0; i < dd.length; i++) krad[i] = Integer.valueOf(dd[i]);

        dd = kc_csv.split(","); if (dd.length==0) return;
        kc = new float[dd.length];
        for (int i = 0; i < dd.length; i++) kc[i] = Float.valueOf(dd[i]);

        dd = maxepoch_csv.split(","); if (dd.length==0) return;
        maxepoch = new int[dd.length];
        for (int i = 0; i < dd.length; i++) maxepoch[i] = Integer.valueOf(dd[i]);

        //******************************************************************
        ImageStack  is_tness;
        ImagePlus   ip_tness=null;

        if (usetness) {

            IJ.log(" -- prefiltering...");
            is_tness = new ImageStack(N, M);
            for (int i = 0; i < P; i++) {
                float[] tt = new float[N * M];
                Arrays.fill(tt, 0f);
                is_tness.addSlice(new FloatProcessor(N, M, tt));
            }

            ip_tness = new ImagePlus("tness", is_tness);
            String[] 	readLn = 	sigmas.trim().split(",");

            for (int i = 0; i < readLn.length; i++) {
                float sig = Float.valueOf(readLn[i].trim()).floatValue();
                TubenessProcessor tp = new TubenessProcessor(sig, false);
                ImagePlus result = tp.generateImage(ip_load);
                IJ.run(result, "Multiply...", "value=" + IJ.d2s(1f/readLn.length,3) + " stack");
                ImageCalculator ic = new ImageCalculator();
                ic.run("Add 32-bit stack", ip_tness, result); // result of the addition is placed in ip_tness
            }

            ip_tness.setCalibration(null);

            if (savemidres) {
                ImagePlus temp = ip_tness.duplicate();
                IJ.run(temp, "8-bit", ""); // convert to 8 bit before saving
                IJ.log("saving... " + midresdir + File.separator + "tness," + sigmas + ".tif");
                IJ.saveAs(temp, "Tiff", midresdir + File.separator + "tness," + sigmas + ".tif");
            }
        }

        // tubeness min-max normalize and store in an array for later and extract locations in a separate array
        float tnessmin = Float.POSITIVE_INFINITY;
        float tnessmax = Float.NEGATIVE_INFINITY;
        tness = new float[SZ];
        int[][] locationXYZ = new int[SZ][3]; // random sampling weighted with the normalized tness as importance function

        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

            float[] slc_float=null;
            byte[] slc_byte=null;

            if (usetness)
                slc_float = (float[]) ip_tness.getStack().getPixels(z);
            else
                slc_byte = (byte[]) ip_load.getStack().getPixels(z);

            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {
                    int ii = (z-1)*(N*M)+y*N+x;

                    if (usetness)
                        tness[ii] = slc_float[y * N + x];
                    else
                        tness[ii] = (float) (slc_byte[y * N + x] & 0xff);

                    locationXYZ[ii][0] = x;
                    locationXYZ[ii][1] = y;
                    locationXYZ[ii][2] = (z-1);

                    if (tness[ii]<tnessmin) tnessmin = tness[ii];
                    if (tness[ii]>tnessmax) tnessmax = tness[ii];

                }
            }
        }

        for (int i = 0; i < SZ; i++) {
            tness[i] = (tnessmax-tnessmin > Float.MIN_VALUE)? ((tness[i]-tnessmin)/(tnessmax-tnessmin)) : 0 ;
        }

        suppmap = new int[SZ];         // suppression map (to avoid re-tracking)
        // go through comma separated parameter values
        for (int i01 = 0; i01 < no.length; i01++) {
            for (int i02 = 0; i02 < ro.length; i02++) {
                for (int i03 = 0; i03 < ni.length; i03++) {
                    for (int i04 = 0; i04 < krad.length; i04++) {
                        for (int i05 = 0; i05 <step.length; i05++) {
                            for (int i06 = 0; i06 < kappa.length; i06++) {
                                for (int i07 = 0; i07 < pS.length; i07++) {
                                    for (int i08 = 0; i08 < pD.length; i08++) {
                                        for (int i09 = 0; i09 < th.length; i09++) {
                                            for (int i10 = 0; i10 < kc.length; i10++) {
                                                for (int i11 = 0; i11 < maxepoch.length; i11++) {

                                                    IJ.log("=================================");
                                                    if (savemidres) {
                                                        Tools.createAndCleanDir(midresdir + File.separator + "g(z|x)");
                                                        Tools.createAndCleanDir(midresdir + File.separator + "suppmap");

                                                        X_swclog = midresdir + File.separator + "Xk.swc";
                                                        X_cnt[0] = 0;
                                                        Tools.cleanfile(X_swclog);
                                                        XP_swclog = midresdir + File.separator + "XPk.swc";
                                                        XP_cnt[0] = 0;
                                                        Tools.cleanfile(XP_swclog);
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

                                                    ImagePlus impool = (usetness)?ip_tness.duplicate():ip_load.duplicate();
                                                    IJ.run(impool, "8-bit", "");

                                                    int threshold = (int) Math.ceil(th[i09] * 255);
                                                    applythreshold(threshold, impool); // will modify values in itnesscopy

                                                    if (savemidres) {
                                                        IJ.saveAs(impool, "Tiff", midresdir + File.separator + "fg,th=" + IJ.d2s(th[i09], 2) + ".tif");
                                                    }

                                                    ArrayList<Integer> locs = new ArrayList<Integer>();     // list of candidate locations for seed points
                                                    ArrayList<Float> locsw = new ArrayList<Float>();       // weights assigned to each location

                                                    for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1
                                                        byte[] slc = (byte[]) impool.getStack().getPixels(z);
                                                        for (int x = 0; x < N; x++) {
                                                            for (int y = 0; y < M; y++) {

                                                                int ii = (z - 1) * (N * M) + y * N + x;

                                                                if ((slc[y * N + x] & 0xff) == 255) {
                                                                    locs.add(ii);
                                                                    locsw.add((float) Math.pow(tness[ii], MultiTT.weight_deg77)); // ((float) Math.exp(tness[ii]));
                                                                }
                                                            }
                                                        }
                                                    }

                                                    if (locs.size() == 0) {
                                                        IJ.log("0 seed candidate locations. no initiation for threshold=" + IJ.d2s(th[i09], 2));
                                                        continue; // try another param configuration
                                                    }

                                                    float locs_count = locs.size();

                                                    if (false && savemidres) {
                                                        // convert before exporting
                                                        ArrayList<int[]> tt = new ArrayList<int[]>(locs.size());
                                                        for (int i = 0; i < locs.size(); i++) {
                                                            int x = locs.get(i) % N;
                                                            int z = locs.get(i) / (N * M);
                                                            int y = locs.get(i) / N - z * M;
                                                            tt.add(new int[]{x, y, z});
                                                        }
                                                        exportlocsxyz(tt, 0.3f, VIOLET, midresdir, "seed_pool");
                                                        tt.clear();
                                                    }

                                                    System.gc();

                                                    IJ.log("-- initialize...");

                                                    MultiTT mtt; // multi-object tracker
                                                    mtt = new MultiTT(P == 1, no[i01], ro[i02], ni[i03], krad[i04], step[i05], kappa[i06], pS[i07], pD[i08], th[i09], kc[i10]);

                                                    if (savemidres) {
                                                        mtt.exporttemplates(midresdir);
                                                        Tools.createAndCleanDir(midresdir + File.separator + "mmodel");
                                                        mtt.mm.getModel(midresdir + File.separator + "mmodel");
                                                    }

                                                    Arrays.fill(suppmap, 0); // reset suppression map, it will fill up as epochs advance

                                                    //******************************************************************
                                                    IJ.log("-- multi-object filtering...");
                                                    long t1 = System.currentTimeMillis();
                                                    int epochcnt = 0;

                                                    while (locs.size() > 0 && epochcnt < maxepoch[i11]) {

                                                        epochcnt++;

                                                        IJ.log("e=" + epochcnt + " [" + maxepoch[i11] + "]");

                                                        iter_count = 0;

                                                        for (int i = locs.size() - 1; i >= 0; i--) {
                                                            if (suppmap[locs.get(i)] > 0) {
                                                                locs.remove(i);
                                                                locsw.remove(i);
                                                            }
                                                        }

                                                        IJ.log(IJ.d2s(locs.size() / 1000f, 1) + "k locations " + IJ.d2s((locs.size()/locs_count)*100f, 1) + "% of the initial pool");

                                                        if (locs.size()==0) {
                                                            IJ.log("locs.size()==0");
                                                            break; // go out of while()
                                                        }

//                                                        ArrayList<int[]> N_o = initlocs(no[i01], locs, locsw); // xyz locations

                                                        ArrayList<int[]> N_o = new ArrayList<int[]>();
                                                        mtt._init(no[i01], TUBE_RADIUS, locs, locsw, N_o, img, N, M, P, tness, suppmap);
                                                        if (N_o.size()==0) {
                                                            IJ.log("initialization stopped, |N_o|=0");
                                                            break; // out of while()
                                                        }

                                                        exportlocsxyz(N_o, 10f, RED, ip_load.getOriginalFileInfo().directory + File.separator, IJ.d2s(N_o.size(),0)+"_seeds_");

                                                        if (savemidres) {

                                                            exportlocsxyz(N_o, 15f, RED, midresdir, "seeds,e=" + IJ.d2s(epochcnt, 0));

                                                            exportXYZW(mtt.Xk, midresdir,       "XWinit,epoch=" + IJ.d2s(epochcnt, 0), MAGENTA);      // phd weights
                                                            exportXYZVxyz(mtt.Xk, midresdir,    "XVinit,epoch=" + IJ.d2s(epochcnt, 0), BLUE);     // directions of particles
                                                            exportNodes(mtt.Y, midresdir,       "Yinit,epoch=" + IJ.d2s(epochcnt, 0)); // estimations

                                                            logval(tnessCsvLog, mtt.cluttertness);
                                                            logval(tnessCsvLog, mtt.Xk);
                                                            logval(zsizeCsvLog, no[i01]);                     // nr. observations mtt.Zk.size()
                                                            logval(phdmassCsvLog, mtt.phdmass);

                                                            if (false) { // it can log the tags... slows down a lot... and takes memory!!!
                                                                String name = "suppmap,_init";
                                                                ImagePlus hpimp = getSuppMap(name);
                                                                IJ.saveAs(hpimp, "Tiff", midresdir + File.separator + "suppmap" + File.separator + name + ".tif");
                                                            }

                                                        }

//                                                        if (true) {IJ.log("blocked iterations"); break;}

                                                        if (mtt.Xk.size() > 0) {

                                                            while (iter_count < maxiter) {

                                                                IJ.log("k=" + iter_count);

                                                                boolean iterok;

                                                                if (iter_count == 0) {
                                                                    iterok = mtt._iter1(N, M, P, tness, suppmap);
//                                                                    iterok = mtt._iter0(N, M, P, tness, suppmap);
                                                                    //iterok = mtt.iter0(N, M, P, tness, suppmap);
                                                                }
                                                                else {
                                                                    iterok = mtt._iter1(N, M, P, tness, suppmap);
                                                                    //iterok = mtt.iter1(N, M, P, tness, suppmap);
                                                                }

                                                                if (savemidres) {

                                                                    Xlog(mtt.Xk, MAGENTA);
                                                                    XPlog(mtt.XPk, WEAK_GREEN);
                                                                    Zlog(mtt.Zk, RED);

                                                                    if (false) { // set if you wish to have the map, takes lot of resources!
                                                                        String name = "suppmap,iter0=" + IJ.d2s(iter_count, 0);
                                                                        ImagePlus hpimp = getSuppMap(name);
                                                                        IJ.saveAs(hpimp, "Tiff", midresdir + File.separator + "suppmap" + File.separator + name + ".tif");
                                                                    }

                                                                    logval(tnessCsvLog, mtt.XPk);
                                                                    logval(zsizeCsvLog, mtt.Zk.size());
                                                                    logval(phdmassCsvLog, mtt.phdmass);

                                                                    if (mtt.g != null) {
                                                                        ImagePlus gimp = new ImagePlus("g(z|x),iter0=" + IJ.d2s(iter_count, 0), new FloatProcessor(mtt.g));
                                                                        IJ.run(gimp, "Rotate 90 Degrees Right", "");
                                                                        IJ.saveAs(gimp, "Tiff", midresdir + File.separator + "g(z|x)" + File.separator + "g(z|x),iter0=" + IJ.d2s(iter_count, 0) + ".tif");
                                                                    }

                                                                }

                                                                if (!iterok) break; // go out of while loop and try new set of seed points

                                                                iter_count++;

                                                            }
                                                        }
                                                        else IJ.log("mtt.Xk.size() == 0");

                                                        // save output before switching to new epoch
                                                        String delindir = imdir + "NDLN.sig.th.no.ro.ni.krad.stp.kapa.ps.pd.kc.e_" + sigmas + "_" + IJ.d2s(th[i09], 2) + "_" + IJ.d2s(no[i01], 0) + "_" + IJ.d2s(ro[i02], 0) + "_" + IJ.d2s(ni[i03], 0) + "_" + IJ.d2s(krad[i04], 0) + "_" + IJ.d2s(step[i05], 0) + "_" + IJ.d2s(kappa[i06], 1) + "_" + IJ.d2s(pS[i07], 2) + "_" + IJ.d2s(pD[i08], 2) + "_" + IJ.d2s(kc[i10], 1) + "_" + IJ.d2s(epochcnt, 0) + "_" + IJ.d2s(new Random().nextInt(Integer.MAX_VALUE),0);
                                                        Tools.createAndCleanDir(delindir);
                                                        remove_double_links(mtt.Y);
                                                        if (!is_biderectinal_linking(mtt.Y)) {
                                                            IJ.log("missing link! fault in bidirectional linking");
                                                            return;
                                                        }

                                                        ArrayList<Node> tree = mtt.bfs1(mtt.Y);
                                                        exportReconstruction(tree, delindir, imnameshort); // export .swc   epochcnt, mtt.Y.size()-1,

                                                    } // while there are locations and epochs have not reached the limit

                                                    long t2 = System.currentTimeMillis();
                                                    IJ.log("done. " + IJ.d2s((t2 - t1) / 1000f, 2) + "s. [maxiter=" + maxiter + ", maxepoch=" + maxepoch[i11] + "]");

                                                    if (false & savemidres) {
                                                        exportDelineation(mtt.Y, midresdir, imnameshort);   // .ndln file
                                                        exportNodes(mtt.Y, midresdir, imnameshort);         // .swc file with isolated nodes
                                                    }

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

    }

    private void remove_double_links(ArrayList<Node> nlist) {
        for (int i = 0; i < nlist.size(); i++) {
            if (nlist.get(i)!=null) {
                Set<Integer> set = new HashSet<Integer>();
                set.addAll(nlist.get(i).nbr);
                nlist.get(i).nbr.clear();
                nlist.get(i).nbr.addAll(set);
            }
        }
    }

    private boolean is_biderectinal_linking(ArrayList<Node> nlist){
        for (int i = 0; i < nlist.size(); i++) {
            if(nlist.get(i)!=null) {
                for (int j = 0; j < nlist.get(i).nbr.size(); j++) {
                    int nbr_idx = nlist.get(i).nbr.get(j);
                    if (Collections.frequency(nlist.get(nbr_idx).nbr, i)!=1) {
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

    private void exportlocsxyz(ArrayList<int[]> locsxyz, float radius, int swctype, String outdir, String swcname) {

        String outfile = outdir + File.separator + swcname + ".swc";

        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            for (int i = 0; i < locsxyz.size(); i++) {
                logWriter.println((i+1) + " " +
                        swctype + " " +
                        IJ.d2s(locsxyz.get(i)[0], 4) + " " +
                        IJ.d2s(locsxyz.get(i)[1], 4) + " " +
                        IJ.d2s(locsxyz.get(i)[2], 4) + " " +
                        IJ.d2s(radius,1) + " " + -1);
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
                float dd = 5; // Xobj.get(i).sig*5f;
                logWriter.println((++swci) + " " + swctype + " " +
                        IJ.d2s(Xobj.get(i).x+ dd *Xobj.get(i).vx,3) + " " +
                        IJ.d2s(Xobj.get(i).y+ dd *Xobj.get(i).vy,3) + " " +
                        IJ.d2s(Xobj.get(i).z+ dd *Xobj.get(i).vz,3) + " " + 0.25 + " " + swcroot);
            }

            logWriter.close();

//            IJ.log("exported : " + outfile);

        } catch (IOException e) {}

    }

//    private void exportXYZSig(ArrayList<X> Xlist, String outdir, String swcname, int type) {
//
//        // exports only locations into swc format for visualization - the rest of the Xk instance is ignored
//        if (outdir==null || swcname==null) return;
//
//        String outfile = outdir + File.separator + swcname + ".swc";
//
//        Tools.cleanfile(outfile);
//
//        try {
//            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));
//
//            logWriter.println("#Xk,Yk,Z,Sig");
//
//            for (int i = 0; i < Xlist.size(); i++) {
//                logWriter.println((i+1) + " " + type + " " +
//                        IJ.d2s(Xlist.get(i).x, 4) + " " +
//                        IJ.d2s(Xlist.get(i).y, 4) + " " +
//                        IJ.d2s(Xlist.get(i).z, 4) + " " +
//                        Xlist.get(i).sig + " " + -1);
//            }
//
//            logWriter.close();
//
//        } catch (IOException e) {}
//
//    }

    private void exportXYZW(ArrayList<X> Xlist, String outdir, String swcname, int type) {

        String outfile = outdir + File.separator + swcname + ".swc";
        Tools.cleanfile(outfile);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

//            logWriter.println("#Xk,Yk,Z,W");

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

    private void exportNodes(ArrayList<Node> nlist, String outdir, String outname) {

        Tools.createDir(outdir);
        String recswc = outdir + File.separator + outname + ".swc";
        Tools.cleanfile(recswc);

        try {

            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(recswc, true)));

            int t = 0;

            for (int i = 0; i < nlist.size(); i++) {
                if (nlist.get(i)!=null) {
                    logWriter.println((++t) + " " + YELLOW + " " +
                            IJ.d2s(nlist.get(i).loc[0], 3) + " " +
                            IJ.d2s(nlist.get(i).loc[1], 3) + " " +
                            IJ.d2s(nlist.get(i).loc[2], 3) + " " +
                            IJ.d2s(nlist.get(i).r, 3) + " " + (-1));
                }
            }

            logWriter.close();

        } catch (IOException e) {}

    }

    private void exportDelineation(ArrayList<Node> nlist, String outdir, String outfile) {

        Tools.createDir(outdir);
        String delinswc1 = outdir + File.separator + outfile + ".ndln";
        Tools.cleanfile(delinswc1);

        try {

            PrintWriter logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(delinswc1, true)));

            int t1 = 0;

            for (int i = 0; i < nlist.size(); i++) {
                if (nlist.get(i)!=null) {

                    logWriter1.print((++t1) + " " + YELLOW + " " + IJ.d2s(nlist.get(i).loc[0], 3) + " " + IJ.d2s(nlist.get(i).loc[1], 3) + " " + IJ.d2s(nlist.get(i).loc[2], 3) + " " + IJ.d2s(nlist.get(i).r, 3) + " ");
                    for (int j = 0; j < nlist.get(i).nbr.size(); j++) {
                        logWriter1.print(IJ.d2s(nlist.get(i).nbr.get(j),0)+" ");
                    }
                    logWriter1.println("");

                }
            }

            logWriter1.close();

        } catch (IOException e) {}

    }

    private void exportReconstruction(ArrayList<Node> nlist, String outdir, String outname) {

        Tools.createDir(outdir);
        String recswc = outdir + File.separator + outname + ".swc";
        Tools.cleanfile(recswc);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(recswc, true)));

            for (int i = 0; i < nlist.size(); i++) {
                if (nlist.get(i)!=null) {

                    Node nn = nlist.get(i);

                    String out =
                            IJ.d2s(i,0)          + " " +
                                    IJ.d2s(nn.type,0)    + " " +
                                    IJ.d2s(nn.loc[0], 3) + " " +
                                    IJ.d2s(nn.loc[1], 3) + " " +
                                    IJ.d2s(nn.loc[2], 3) + " " +
                                    IJ.d2s(nn.r,      3) + " " +
                                    ((nn.nbr.size()==0)?"-1":IJ.d2s(nn.nbr.get(0),0));

                    logWriter.println(out);

                    if (nn.nbr.size()>1)
                        IJ.log("*** ERROR in tree export "+i);
                }
            }

            logWriter.close();

        } catch (IOException e) {}

    }

    private void exportReconstruction(int numepochs, int numnodes, ArrayList<Node> nlist, String outdir, String outname) {//

        outdir += IJ.d2s(numepochs,0)+"_"+IJ.d2s(numnodes,0);
        Tools.createDir(outdir);
        String recswc = outdir + File.separator + outname + ".swc";
        Tools.cleanfile(recswc);

        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(recswc, true)));

            for (int i = 0; i < nlist.size(); i++) {
                if (nlist.get(i)!=null) {

                    Node nn = nlist.get(i);

                    String out =
                            IJ.d2s(i,0)          + " " +
                            IJ.d2s(nn.type,0)    + " " +
                            IJ.d2s(nn.loc[0], 3) + " " +
                            IJ.d2s(nn.loc[1], 3) + " " +
                            IJ.d2s(nn.loc[2], 3) + " " +
                            IJ.d2s(nn.r,      3) + " " +
                                    ((nn.nbr.size()==0)?"-1":IJ.d2s(nn.nbr.get(0),0));

                    logWriter.println(out);

                    if (nn.nbr.size()>1)
                        IJ.log("*** ERROR in tree export "+i);
                }
            }

            logWriter.close();

        } catch (IOException e) {}

    }

    public ImagePlus getSuppMap(String title){

        ImageStack outis = new ImageStack(N, M);

        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

            float[] slc = new float[N*M];

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

    private int[] gethist(ImagePlus imp, boolean usestack) {

        ImageProcessor ip = imp.getProcessor();

        ImageStatistics stats = null;
        if (usestack) {
            stats = new StackStatistics(imp);
        }
        else {
            if (!(ip instanceof ByteProcessor)) {
                ip.resetMinAndMax();
                imp.updateAndDraw();
            }
            stats = ImageStatistics.getStatistics(ip, Measurements.AREA+ Measurements.MIN_MAX+Measurements.MODE, null);
        }

        int maxCount2 = 0;
        int[] histogram = new int[stats.nBins];

        for (int i = 0; i < stats.nBins; i++)
            histogram[i] = stats.histogram[i];

        for (int i = 0; i < stats.nBins; i++)
            if ((histogram[i] > maxCount2) && (i != stats.mode))
                maxCount2 = histogram[i];

        int hmax = stats.maxCount;

        if ((hmax>(maxCount2 * 2)) && (maxCount2 != 0)) {
            hmax = (int)(maxCount2 * 1.5);
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
            for (int j = 0; j < w*h; j++) {
                if ((lay[j] & 0xff) >= thval)   lay[j] = (byte) 255;
                else                            lay[j] = (byte) 0;
            }
        }

    }

    private ArrayList<Integer> importsamp(ArrayList<Double> lcws, int n) {
        // systematic resampling, Beyond Kalman Filtering, Ristic et al.
        double totalmass = lcws.get(lcws.size()-1);
        double u1 = (totalmass/(float)n) * new Random().nextDouble();

        ArrayList<Integer> out = new ArrayList<Integer>(n);
        out.clear();
        int i = 0;
        for (int j = 0; j < n; j++) {
            double uj = u1 + j*(totalmass/(float)n);
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
                    if (Math.pow(x,2)/Math.pow(rxy,2) + Math.pow(y,2)/Math.pow(rxy,2) +
                            ((rz>0)?Math.pow(z,2)/Math.pow(rz,2):0)
                            <= 1) {
                        px.add(x); py.add(y); pz.add(z);
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

//    private ArrayList<int[]> initlocs(int n, ArrayList<Integer> l, ArrayList<Float> w){}

}
