package com.braincadet.ndelin.multi;

import com.braincadet.ndelin.fun.Bessel;
import com.braincadet.ndelin.fun.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.ArrayList;

/**
 * Stepper class has some commonly used offsets and weights and the directions assigned to the offsets
 * there is a set of weights used for importance sampling for object particle prediction
 * also suppression weight function (included in importance sampling)
 * Created by miroslav on 11/26/15.
 */
public class Stepper {

    boolean is2d;
    int R2;                              // outer radius
    int R1;                              // inner radius
    public static int ndirs2d = 36;      //
    public static int ndirs3d = 60;      //
    float kappa;
    float sig;

    public int[][]    p;        // offs2 * 3
    public float[]    pcws;     // storage for the cummulative sum, does not have to contain all values
    public double[]   d;        // offs2
    public double[][] u;        // offs2 * 3

    public int [][]   psup;    // offs1 * 3
    public double[]   wsup;    // offs1
    
    public double[][] w;        // ndir * offs
    public double[][] v;        // ndir * 3

    public int cloudxyz[][]; // could of predefined number of points filling up the sphere

    public Stepper(int step, boolean is2d, float kappa, int ro) {

        this.R2 = 2*step;
        this.R1 = step;// (step/2<1)?1:(step/2);
        this.is2d = is2d;
        this.kappa = kappa;
        this.sig = step/2f; // by convention

        ArrayList<int[]> off2 = new ArrayList<int[]>(); // R2
        ArrayList<int[]> off1 = new ArrayList<int[]>(); // R1
        
        for (int x = -R2; x <= R2; x++) {
            for (int y = -R2; y <= R2; y++) {

                if (is2d) {

                    int d2 = x*x+y*y;

                    if (d2>0 && d2<=R2*R2)
                        off2.add(new int[]{x, y, 0});

                    if (d2>=0 && d2<=R1*R1)
                        off1.add(new int[]{x, y, 0});

                }
                else {

                    for (int z = -R2; z <= R2; z++) {

                        int d2 = x*x+y*y+z*z;

                        if (d2>0 && d2<= R2*R2)
                            off2.add(new int[]{x, y, z});

                        if (d2>=0 && d2<=R1*R1)
                            off1.add(new int[]{x, y, z});

                    }
                }
            }
        }

        p = new int[off2.size()][3];
        pcws = new float[off2.size()];
        d = new double[off2.size()];
        u = new double[off2.size()][3];

        for (int i = 0; i < off2.size(); i++) {

            p[i][0] = off2.get(i)[0];
            p[i][1] = off2.get(i)[1];
            p[i][2] = off2.get(i)[2];

            d[i] = Math.sqrt(Math.pow(p[i][0],2)+Math.pow(p[i][1],2)+Math.pow(p[i][2],2));

            u[i][0] = p[i][0]/d[i];
            u[i][1] = p[i][1]/d[i];
            u[i][2] = p[i][2]/d[i];
        }

        psup = new int[off1.size()][3];
        wsup = new double[off1.size()];

        for (int i = 0; i < off1.size(); i++) {

            psup[i][0] = off1.get(i)[0];
            psup[i][1] = off1.get(i)[1];
            psup[i][2] = off1.get(i)[2];

            wsup[i] = Math.pow(  Math.pow(psup[i][0],2)+Math.pow(psup[i][1],2)+Math.pow(psup[i][2],2)  , 3f/2);

        }

        v = Tools.calcdirs360(this.is2d, (this.is2d ? ndirs2d : ndirs3d)); // form the directions

        double rad, circ, val, dotp; // radial, polar distance

        w = new double[v.length][off2.size()];

        for (int i = 0; i < v.length; i++) { // v.length == number of directions

            float wmin=Float.POSITIVE_INFINITY, wmax=Float.NEGATIVE_INFINITY;
            float wsum = 0;

            for (int j = 0; j < p.length; j++) {

                rad = Math.exp(-Math.pow(d[j]-R1, 2)/(2*Math.pow(sig,2)));

                dotp = v[i][0]*u[j][0] + v[i][1]*u[j][1] + v[i][2]*u[j][2];
                dotp = (dotp>1)? 1 : (dotp<-1)? -1 : dotp;
                circ = Math.exp(kappa * dotp) / (2.0*3.14* Bessel.I0(kappa)); // von misses distribution

                val = circ * rad;
                w[i][j] = (float) val;

                wsum += val;

                if (val<wmin) wmin = (float) val;
                if (val>wmax) wmax = (float) val;

            }

//            IJ.log("v["+i+"], wmin="+wmin+", wmax="+wmax + " " + Arrays.toString(v[i]));

            for (int j = 0; j < p.length; j++)
                w[i][j] = w[i][j]/wsum;
//                w[i][j] = (w[i][j] - wmin) / (wmax - wmin);

        }

        //***************************
//        int ro = step/2;
//        ro = (ro<1)? 1 : ro;
        int T = (int) Math.ceil((is2d)?Math.pow(ro, 1f/2):Math.pow(ro,1f/3))/2; // +/- T

        ArrayList<int[]> offs = new ArrayList<int[]>();
        ArrayList<Float> dist2 = new ArrayList<Float>();
        for (int dx = -T; dx <= T; dx++) {
            for (int dy = -T; dy <= T; dy++) {
                if (is2d) {
                    offs.add(new int[]{dx, dy, 0});
                    dist2.add((float) (dx*dx+dy*dy+0*0));
                }
                else {
                    for (int dz = -T; dz <= T; dz++) {
                        offs.add(new int[]{dx, dy, dz});
                        dist2.add((float) (dx*dx+dy*dy+dz*dz));
                    }
                }
            }
        }

        int[] idxs = Tools.asc(dist2); // sort

        cloudxyz = new int[ro][3];

        for (int i = 0; i < ro; i++) {

            cloudxyz[i][0] = offs.get(idxs[i])[0];
            cloudxyz[i][1] = offs.get(idxs[i])[1];
            cloudxyz[i][2] = offs.get(idxs[i])[2];

        }

    }

    public int getdirection(float vx, float vy, float vz) {

        int idx = -1;
        double maxdotp = Double.NEGATIVE_INFINITY;
        double currdotp;
        for (int i = 0; i < v.length; i++) {

            currdotp = vx* v[i][0] + vy* v[i][1] + vz* v[i][2];
            if (currdotp>maxdotp) {
                maxdotp = currdotp;
                idx = i;
            }

        }

        return idx;

    }

    public void getModel(String outdir) {

        // direction sampling
        String name = outdir + File.separator + "dirs.swc";
        Tools.cleanfile(name);
        try {
            PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(name, true)));
            int swci = 0;
            logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(0,2) + " " + IJ.d2s(0,2) + " " + IJ.d2s(0,2) + " " + 0.05 + " " + (-1));
            int swcroot = swci;

            for (int j = 0; j < v.length; j++) {
                logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(R2 * v[j][0],3) + " " + IJ.d2s(R2 * v[j][1],3) + " " + IJ.d2s(R2 * v[j][2],3) + " " + 0.05 + " " + swcroot);
            }

            logWriter.close();

        } catch (IOException e) {}

        // export swc and tif of the wieghts corresponding to each direction
        for (int i = 0; i < w.length; i++) {

            // save swc
            name = outdir + File.separator + "dir" + IJ.d2s(i,0) + ",v=" + IJ.d2s(v[i][0],2) + "," + IJ.d2s(v[i][1],2) + "," + IJ.d2s(v[i][2],2)+".swc";
            Tools.cleanfile(name);
            try {
                PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(name, true)));
                int swci = 0;
                for (int j = 0; j < w[i].length; j++) {
                    logWriter.println((++swci) + " " + 2 + " " + IJ.d2s(p[j][0],3) + " " + IJ.d2s(p[j][1],3) + " " + IJ.d2s(p[j][2],3) + " " + IJ.d2s(0.5f* w[i][j],3) + " " + -1);
                }
                logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(0,3) + " " + IJ.d2s(0,3) + " " + IJ.d2s(0,3) + " " + 1 + " " + -1);
                logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(0,3) + " " + IJ.d2s(0,3) + " " + IJ.d2s(0,3) + " " + 0.25 + " " + -1);
                int swcroot = swci;
                logWriter.println((++swci) + " " + 1 + " " + IJ.d2s(R2 * v[i][0],3) + " " + IJ.d2s(R2 * v[i][1],3) + " " + IJ.d2s(R2 * v[i][2],3) + " " + 0.25 + " " + swcroot);

                logWriter.close();

            } catch (IOException e) {}

            // save image/image stack for offsets and directions for particular direction
            ImageStack isout = new ImageStack(2*R2+1, 2*R2+1);

            float[][] t = new float[(is2d)?1:(2*R2+1)][(2*R2+1)*(2*R2+1)];

            for (int j = 0; j < w[i].length; j++) {

                int x = p[j][0] + R2;
                int y = p[j][1] + R2;
                int z = p[j][2] +((is2d)?0: R2);

                t[z][y*(2* R2 +1)+x] = (float) w[i][j];

            }

            for (int j = 0; j < t.length; j++) {

                FloatProcessor fpout = new FloatProcessor(2* R2 +1, 2* R2 +1, t[j]);
                isout.addSlice(fpout);

            }

            name = outdir + File.separator + "dir" + IJ.d2s(i,0) + ",v=" + IJ.d2s(v[i][0],2) + "," + IJ.d2s(v[i][1],2) + "," + IJ.d2s(v[i][2],2)+".tif";

            ImagePlus ipout = new ImagePlus(name, isout);
//            IJ.run(ipout, "8-bit", "");
            IJ.saveAs(ipout, "Tiff", name);

        }

        // save image/image stack for offsets and directions for particular direction
        ImageStack isout = new ImageStack(2*R1+1, 2*R1+1);

        float[][] t = new float[(is2d)?1:(2*R1+1)][(2*R1+1)*(2*R1+1)];

        for (int j = 0; j < psup.length; j++) {

            int x = psup[j][0] + R1;
            int y = psup[j][1] + R1;
            int z = psup[j][2] +((is2d)?0: R1);

            t[z][y*(2*R1+1)+x] = (float) wsup[j];

        }

        for (int j = 0; j < t.length; j++) {

            FloatProcessor fpout = new FloatProcessor(2*R1+1, 2*R1+1, t[j]);
            isout.addSlice(fpout);

        }

        name = outdir + File.separator + "supp.tif";
        ImagePlus ipout = new ImagePlus(name, isout);
//            IJ.run(ipout, "8-bit", "");
        IJ.saveAs(ipout, "Tiff", name);

    }

}
