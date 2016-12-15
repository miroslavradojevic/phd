package com.braincadet.phd.single;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.File;

/*
Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
Created by miroslav on 9/21/15.
*/

public class STracker implements PlugIn, MouseListener, MouseMotionListener {

    ImageCanvas cnv;
    int     sc      = 10;       // scale (~neurite diameter)
    float   angdeg  = 40;       // angular standard deviation in degrees
    int     ni      = 20;       // max nr iterations
    int     np      = 50;       // number of particles

    float[] img;
    int N, M, P;
    String outputdir;

    // object tracker classes
    SingleTT stt;

    public void mouseClicked(MouseEvent mouseEvent) {

        int x = cnv.offScreenX(mouseEvent.getX());
        int y = cnv.offScreenY(mouseEvent.getY());
        int z = cnv.getImage().getZ()-1;

        IJ.log(x + ", " + y + ", " + z + " I=" + (img[z*(N*M)+y*N+x]));

        float[] vxyz = stt.locdirection(img, N, M, P, x, y, z, outputdir);

        stt.track(img, N, M, P, x, y, z, vxyz[0], vxyz[1], vxyz[2], sc, angdeg*(3.14f/180f), ni, np, outputdir);

    }

    public void mousePressed(MouseEvent mouseEvent) {}
    public void mouseReleased(MouseEvent mouseEvent) {}
    public void mouseEntered(MouseEvent mouseEvent) {}
    public void mouseExited(MouseEvent mouseEvent) {}
    public void mouseDragged(MouseEvent mouseEvent) {}
    public void mouseMoved(MouseEvent mouseEvent) {}

    public void run(String s) {

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory();
        String image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus ip_load = new ImagePlus(image_path);
        if(ip_load==null) return;
        if(ip_load.getType()!=ImagePlus.GRAY8) {IJ.log("image needs to be gray8"); return;}

        N = ip_load.getWidth();
        M = ip_load.getHeight();
        P = ip_load.getStack().getSize();

        ip_load.setCalibration(null);

        // read image array
        img = new float[N*M*P];
        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

            byte[] slc = (byte[]) ip_load.getStack().getPixels(z);

            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {
                    img[(z-1)*(N*M)+y*N+x] = slc[y*N+x] & 0xff;
                }
            }

        }

        // load tracking parameters
        GenericDialog gd = new GenericDialog("SMC-PF");
        gd.addNumericField("SCALE",          sc,     0, 5, "pix");
        gd.addNumericField("ANG",            angdeg, 0, 5, "deg");
        gd.addNumericField("ITER",           ni,     0, 5, "");
        gd.addNumericField("NPARTICLES", 	 np,     0, 5, "");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        sc      = (int) gd.getNextNumber();
        angdeg  = (float) gd.getNextNumber();
        ni	    = (int) gd.getNextNumber();
        np      = (int) gd.getNextNumber();

        // initialize tracker
        stt = new SingleTT(sc, P==1);
        outputdir = ip_load.getOriginalFileInfo().directory+ip_load.getTitle()+"_out";
        File f1 = new File(outputdir);
        if (!f1.exists()) {
            f1.mkdirs();
        }
        else { // clean it, reset with each exec
            File[] files = f1.listFiles();
            for (File file : files)
            {
                if (!file.delete())
                    System.out.println("Failed to delete " + file);
            }
        }

        IJ.log("exporttemplates...");
        stt.exporttemplates(outputdir);
        IJ.log("demo...");
        stt.demo(np, 5, angdeg*(3.14f/180f), outputdir);
        IJ.log("done.");

        ip_load.show();
        cnv = ip_load.getCanvas();
        cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);
        IJ.setTool("hand");

    }

}
