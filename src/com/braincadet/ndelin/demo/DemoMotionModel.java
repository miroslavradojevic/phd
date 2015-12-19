package com.braincadet.ndelin.demo;

import com.braincadet.ndelin.fun.Tools;
import com.braincadet.ndelin.multi.Stepper;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.io.File;
import java.util.Random;

/**
 * Created by miroslav on 11/26/15.
 */
public class DemoMotionModel implements PlugIn {

    public void run(String s) {

        IJ.log("DemoMotionModel...");

        int step = 5;
        boolean is2d = false;
        float kappa = 1;
        String outdir = System.getProperty("user.home") + File.separator + "mmodel";

        // load tracking parameters
        GenericDialog gd = new GenericDialog("MMODEL");
        gd.addNumericField("step",      step,   0, 5, "pix");
        gd.addNumericField("kappa",     kappa,  0, 5, "");
        gd.addCheckbox("2d",            is2d);
        gd.addStringField("outdir",     outdir, 50);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        step        = (int) gd.getNextNumber();
        kappa       = (float) gd.getNextNumber();
        is2d	    = gd.getNextBoolean();
        outdir      = gd.getNextString();

        Tools.createAndCleanDir(outdir);

        // generate random direction
        float vx = new Random().nextFloat() * 2f - 1f;
        float vy = new Random().nextFloat() * 2f - 1f;
        float vz = new Random().nextFloat() * 2f - 1f;
        vz = (is2d)? 0 : vz;
        float vn = (float) Math.sqrt(vx * vx + vy * vy + ((is2d)?0:(vz * vz)));
        vx /= vn;   vy /= vn;   vz = (is2d)? 0 : vz/vn;

        // usage of the Stepper class
        Stepper mm = new Stepper(step, is2d, kappa, 10);
        mm.getModel(outdir);
        int idx = mm.getdirection(vx, vy, vz);

        IJ.log(IJ.d2s(vx,2)+","+IJ.d2s(vy,2)+","+IJ.d2s(vz,2) + " -- I="
                + idx + "("+((is2d)?mm.ndirs2d:mm.ndirs3d)+"):" + IJ.d2s(mm.v[idx][0],2)+","+IJ.d2s(mm.v[idx][1],2)+","+IJ.d2s(mm.v[idx][2],2));
        IJ.log("done.");

    }
}
