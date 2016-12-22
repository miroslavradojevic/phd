package com.braincadet.phd;

import com.braincadet.phd.Bessel;
import com.braincadet.phd.Tools;
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
//    float sig;

    public int        sz;
    public int[][]    p;        // offs2 * 3
    public float[]    pcws;     // storage for the cummulative sum, does not have to contain all values
    public double[]   d;        // offs2
    public double[][] u;        // offs2 * 3
    public double[][] w;        // ndir * offs
    // perhaps add w_cws to be complete, for now, it's not necessary
    public double[]   w0;       // offs
    public double[][] v;        // ndir * 3

    public int cloudxyz[][];    // could of predefined number of points filling up the sphere
    public int Nxyz[][][];        //

    public int        _sz;      // count
    public int[][]    _p;       // xyz offsets
    public int[][]    _p0;      // xyz offsets
    public float[]    _pcwsX;    // cummulative weight sum for sampling X predicted particles
    public float[]    _pcwsZ;    // cummulative weight sum for sampling
    public float[]    _pcws0;   // cummulative weight sum (center included)
    public double[]   _d;       //
    public double[][] _u;       //
    public double[][] _w;       //
    public double[][] _w_cws;   //
    public double[]   _w0;      //

    public Stepper(int step, boolean is2d, float kappa, int ro) {

        // ro will initialize T, cloudxyz and Nxyz

        this.R2 = 2*step;
        this.R1 = step;// (step/2<1)?1:(step/2);
        this.is2d = is2d;
        this.kappa = kappa;

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

        sz = off2.size();
        p = new int[sz][3];
        pcws = new float[sz];
        d = new double[sz];
        u = new double[sz][3];

        for (int i = 0; i < sz; i++) {

            p[i][0] = off2.get(i)[0];
            p[i][1] = off2.get(i)[1];
            p[i][2] = off2.get(i)[2];

            d[i] = Math.sqrt(Math.pow(p[i][0],2)+Math.pow(p[i][1],2)+Math.pow(p[i][2],2));

            u[i][0] = p[i][0]/d[i];
            u[i][1] = p[i][1]/d[i];
            u[i][2] = p[i][2]/d[i];

//            IJ.log("u["+i+"]=("+u[i][0]+","+u[i][1]+","+u[i][2]+")");
        }

        v = Tools.calcdirs360(this.is2d, (this.is2d ? ndirs2d : ndirs3d)); // form the directions

        double rad, circ, val, dotp; // radial, polar distance

        w  = new double[v.length][sz];
        double w_min = Float.POSITIVE_INFINITY;
        double w_max = Float.NEGATIVE_INFINITY;

        for (int i = 0; i < v.length; i++) { // v.length == number of directions



            for (int j = 0; j < sz; j++) {

                rad = Math.exp(-Math.pow(d[j]-R1, 2)/(2*Math.pow(R1/3f,2)));

                dotp = v[i][0]*u[j][0] + v[i][1]*u[j][1] + v[i][2]*u[j][2];
                dotp = (dotp>1)? 1 : (dotp<-1)? -1 : dotp;
                circ = Math.exp(kappa * dotp) / (2.0*3.14* Bessel.I0(kappa)); // von misses distribution

                val = (dotp>0)? circ * rad : 0 ;
                w[i][j] = (float) val;

                if (w[i][j]<w_min) w_min = w[i][j];
                if (w[i][j]>w_max) w_max = w[i][j];

            }

            float wsum = 0;

            for (int j = 0; j < sz; j++) {

                w[i][j] = ((w[i][j]-w_min)/(w_max-w_min)); // z normalized
                w[i][j] = (w[i][j]<0.2)? 0 : w[i][j]; // crop so that there is no backward prediction
                wsum += w[i][j];

            }

            for (int j = 0; j < p.length; j++) {
                w[i][j] = w[i][j] / wsum; // normalize into probability distribution
            }

//            String tt = "";
//            for (int j = 0; j < p.length; j++) {
//                tt += IJ.d2s(((w[i][j]-w_min)/(w_max-w_min)), 4) +" | ";
//            }
//            IJ.log(tt);

        }

        w0 = new double[off2.size()];  // this one has no directions embedded

        float w0sum = 0;

        for (int j = 0; j < p.length; j++) {
            w0[j] = (float) (d[j]/R2);
            w0sum += w0[j];

        }

        for (int j = 0; j < p.length; j++)
            w0[j] = w0[j]/w0sum;

        //******************************************************
        if (ro>0) {

            int T = (int) Math.ceil((is2d)?Math.pow(ro, 1f/2):Math.pow(ro,1f/3))/2; // +/- T  T - radius to accommodate up to ro points in sphere

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
                            dist2.add((float) (dx*dx+dy*dy+(1.5*dz)*(1.5*dz)));
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

            //******************************************************
            Nxyz = new int[ro+1][][];
            Nxyz[0] = null;

            for (int i = 1; i <= ro; i++) {

                Nxyz[i] = new int[i][3];

                for (int k = 0; k < i; k++) {
                    Nxyz[i][k][0] = offs.get(idxs[k])[0];
                    Nxyz[i][k][1] = offs.get(idxs[k])[1];
                    Nxyz[i][k][2] = offs.get(idxs[k])[2];
                }

            }

        }
        else {
            cloudxyz = null;
            Nxyz = null;
        }

        //******************************************************
        ArrayList<int[]> _off = new ArrayList<int[]>();

        for (int x = -R1; x <= R1; x++) {
            for (int y = -R1; y <= R1; y++) {

                if (is2d) {

                    int d2 = x*x+y*y;

                    if (d2>0 && d2<=R1*R1)
                        _off.add(new int[]{x, y, 0});

                }
                else {

                    for (int z = -R1; z <= R1; z++) {

                        int d2 = x*x+y*y+z*z;

                        if (d2>0 && d2<=R1*R1)
                            _off.add(new int[]{x, y, z});

                    }
                }
            }
        }

        _sz = _off.size();

        _p      = new int[_sz][3];
        _p0     = new int[_sz+1][3];
        _pcwsX   = new float[_sz];
        _pcwsZ   = new float[_sz];
        _pcws0  = new float[_sz+1];
        _d      = new double[_sz];
        _u      = new double[_sz][3];

        for (int i = 0; i < _sz; i++) {

            _p0[i][0] = _p[i][0] = _off.get(i)[0];
            _p0[i][1] = _p[i][1] = _off.get(i)[1];
            _p0[i][2] = _p[i][2] = _off.get(i)[2];

            _d[i] = Math.sqrt(Math.pow(_p[i][0],2)+Math.pow(_p[i][1],2)+Math.pow(_p[i][2],2));

            _u[i][0] = _p[i][0]/_d[i];
            _u[i][1] = _p[i][1]/_d[i];
            _u[i][2] = _p[i][2]/_d[i];

        }

        _p0[_sz][0] = _p0[_sz][1] = _p0[_sz][2] = 0;

        _w  = new double[v.length][_sz];
        _w_cws = new double[v.length][_sz]; // cummulative weight sum (used for sampling)

        for (int i = 0; i < v.length; i++) { // v.length == number of directions

            float wsum = 0;

            for (int j = 0; j < _sz; j++) {

                dotp = v[i][0]*_u[j][0] + v[i][1]*_u[j][1] + v[i][2]*_u[j][2];
                dotp = (dotp>1)? 1 : (dotp<-1)? -1 : dotp;
                circ = Math.exp(kappa * dotp) / (2.0*3.14* Bessel.I0(kappa)); // von misses distribution

                val = circ * 1f; // this set had radial component cancelled
                _w[i][j] = (float) val;

                wsum += val;

            }

            for (int j = 0; j < _sz; j++) {
                _w[i][j] = _w[i][j] / wsum;
                _w_cws[i][j] = (j==0)? _w[i][j] : (_w[i][j]+_w_cws[i][j-1]);
            }

//            IJ.log("just check " + _w_cws[i][_sz-1]);

        }

        _w0 = new double[_sz];  // this one has no directions embedded

        float _w0sum = 0;

        for (int j = 0; j < _sz; j++) {
            _w0[j] = 1f; // (float) (d[j]/R2);
            _w0sum += _w0[j];

        }

        for (int j = 0; j < _sz; j++)
            _w0[j] = _w0[j]/_w0sum;

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

        if (idx==-1) {
            IJ.log("---");
            IJ.log("vx,vy,vz: " + vx + "," + vy + "," + vz + "   ");
            IJ.log("---");
            String ss = "";
            for (int i = 0; i < v.length; i++) ss+="[" + v[i][0] + "," + v[i][1] + "," + v[i][2] + "] ";
            IJ.log(ss);
            IJ.log("---");
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

        // export w for each direction
        for (int i = 0; i < v.length; i++) {

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

        // export _w for each direction
        for (int i = 0; i < v.length; i++) {

            ImageStack isout = new ImageStack(2*R1+1, 2*R1+1);

            float[][] t = new float[(is2d)?1:(2*R1+1)][(2*R1+1)*(2*R1+1)];

            for (int j = 0; j < _w[i].length; j++) {

                int x = _p[j][0] + R1;
                int y = _p[j][1] + R1;
                int z = _p[j][2] +((is2d)?0: R1);

                t[z][y*(2* R1 +1)+x] = (float) _w[i][j];

            }

            for (int j = 0; j < t.length; j++) {

                FloatProcessor fpout = new FloatProcessor(2*R1+1, 2*R1+1, t[j]);
                isout.addSlice(fpout);

            }

            name = outdir + File.separator + "_w" + IJ.d2s(i,0) + ",v=" + IJ.d2s(v[i][0],2) + "," + IJ.d2s(v[i][1],2) + "," + IJ.d2s(v[i][2],2)+".tif";

            ImagePlus ipout = new ImagePlus(name, isout);
//            IJ.run(ipout, "8-bit", "");
            IJ.saveAs(ipout, "Tiff", name);

        }

    }

}
