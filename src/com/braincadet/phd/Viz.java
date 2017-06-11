package com.braincadet.phd;

import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;

import java.io.File;

public class Viz implements PlugIn {

    public void run(String s) {

        String p1;
        int p2;

        OpenDialog.setDefaultDirectory(System.getProperty("user.home"));
        OpenDialog dc = new OpenDialog("Select image");
        String image_path = dc.getPath();
        if (image_path == null) return;

        ImagePlus ip = new ImagePlus(image_path); // read selected image stack

        if (Macro.getOptions() == null) {
            //
            GenericDialog gd = new GenericDialog("VIZ");
            gd.addStringField("param1", "dummy", 10);
            gd.addNumericField("param2",  99, 0, 5, "");

            gd.showDialog();
            if (gd.wasCanceled()) return;

            p1 = gd.getNextString();
            p2 = (int) gd.getNextNumber();

        }
        else {
            p1 = Macro.getValue(Macro.getOptions(), "param1", "dummy");
            p2 = Integer.valueOf(Macro.getValue(Macro.getOptions(), "param2", String.valueOf(99)));
        }

        // max along z
        ZProjector zp = new ZProjector(ip);
        zp.setMethod(ZProjector.MAX_METHOD);
        zp.doProjection();
        ImagePlus ip1 = zp.getProjection();
        IJ.run(ip1, "Invert", "");

        createDir(ip.getOriginalFileInfo().directory + "viz");

        IJ.saveAs(ip1, "TIF", ip.getOriginalFileInfo().directory + "viz"+ File.separator+ip1.getTitle());

    }

    private static void createDir(String dirpath) {
        // create directory without cleaning it up
        File f1 = new File(dirpath);
        if (!f1.exists()) {f1.mkdirs();}
    }

}
