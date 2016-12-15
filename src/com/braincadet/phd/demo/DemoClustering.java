package com.braincadet.ndelin.demo;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Vector;

/**
 * Created by miroslav on 10/26/15.
 */
public class DemoClustering implements PlugIn {

    Vector<float[]> disks;
    static int N = 30;

    public static void main (String[] args) {}

    private static Color getRandomColor(){
        int choose = new Random().nextInt(10);
        switch (choose) {
            case 0: return Color.RED;
            case 1: return Color.GREEN;
            case 2: return Color.WHITE;
            case 3: return Color.ORANGE;
            case 4: return Color.PINK;
            case 5: return Color.YELLOW;
            case 6: return Color.GRAY;
            case 7: return Color.MAGENTA;
            case 8: return Color.CYAN;
            case 9: return Color.BLUE;
            default:return Color.BLACK;
        }

    }

    /*
        problem 1: group centroids with different radiuses (disks) to the same cluster if they overlap
        this way each component maintains it's own radius (useful when clustering regions of different disk sizes)
     */
    public static int[] clustering(Vector<float[]> disks) //[x, y, r]  // Vector<float[]>
    {

        int[] labels = new int[disks.size()];
        for (int i = 0; i < labels.length; i++) labels[i] = i;

        IJ.log("INIT. LABELS:");
        for (int i = 0; i < labels.length; i++)
            System.out.print(labels[i]+" ");
        IJ.log("");

        for (int i = 0; i < disks.size(); i++) {

            // one versus the rest
            for (int j = 0; j < disks.size(); j++) {

                if (i!=j) {

                    double dst2 	= Math.pow(disks.get(i)[0]-disks.get(j)[0], 2) + Math.pow(disks.get(i)[1]-disks.get(j)[1], 2);
                    double rd2 		= Math.pow(disks.get(i)[2]+disks.get(j)[2], 2);
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

        IJ.log("OUT LABELS:");
        for (int ii = 0; ii < labels.length; ii++)
            System.out.print(labels[ii]+" ");
        IJ.log("");

        return labels; // cluster labels for each disc

    }

    /*
        problem 2: group centroids with fixed radius used as threshold distance
        (useful when clustering regions with the same disk sizes)
     */
    public static int[] clustering(int[] idxs, int[][] dists, int threshold_dists)
    {
        // indxs represent indexes of values that need to be clustered,
        // dists are the distances,
        // threshold_dists is the distance limit
        // output is list of unique labels
        int[] labels = new int[idxs.length];
        for (int i = 0; i < labels.length; i++) labels[i] = i;

        IJ.log("INIT. LABELS:");
        IJ.log(Arrays.toString(labels));

        for (int i = 0; i < idxs.length; i++) {

            // one versus the rest
            for (int j = 0; j < idxs.length; j++) {

                //
                if (i != j) {

                    if (dists[i][j]<=threshold_dists) {

                        if (labels[j] != labels[i]) {

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

        IJ.log("OUT LABELS:");
//		for (int ii = 0; ii < labels.length; ii++)
//			System.out.print(labels[ii]+" ");
        IJ.log(Arrays.toString(labels));

        return labels;

    }

    // use output of clustering to give out the final cluster centroids
    public static ArrayList<float[]> extracting(int[] labels, int[] vals)
    { // int[] idxs,

        boolean[] checked = new boolean[labels.length];
        ArrayList<float[]> out = new ArrayList<float[]>();

        for (int i = 0; i < labels.length; i++) {
            if (!checked[i]) {

                float centroid = vals[ i ]; // idxs[i]
                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < labels.length; j++) {
                    if (!checked[j]) {
                        if (labels[j]==labels[i]) {

                            centroid += vals[ j ]; // idxs[j]
                            count++;
                            checked[j] = true;

                        }
                    }
                }

                out.add(new float[]{centroid/count, count});

            }
        }

        return out;

    }

    public void run(String s) {

        float minx = 20,    miny = 20;
        float maxx = 100,   maxy = 100;
        float minr = 2,     maxr = 5;

        int H = 120, W = 120;

        // table with colors
        Color[] 		c = new Color[N];
        for (int i = 0; i < N; i++) c[i] = getRandomColor();

        float[][] p = new float[N][];

        // fille them up
        for (int i = 0; i < N; i++) {
            float x = minx + new Random().nextFloat() * (maxx-minx);
            float y = miny + new Random().nextFloat() * (maxy-miny);
            float r = minr + new Random().nextFloat() * (maxr-minr);
            p[i] = new float[]{x, y, r};
        }

        IJ.log("#########################");
        IJ.log("#### TEST 1 ####");
        IJ.log("#########################");

        // input disks 2d
        Vector<float[]> disks = new Vector<float[]>(N);

        for (int i = 0; i < N; i++) {
            disks.add(p[i]);
        }

        IJ.log("\ncreated " + disks.size() + " points:\n");
        int[] lab = clustering(disks);
        IJ.log(Arrays.toString(lab));

        ByteProcessor ip = new ByteProcessor(W, H);
        ImagePlus im = new ImagePlus("clustering2d", ip);

        Overlay ov = new Overlay();
        ov.drawNames(true);

        for (int i = 0; i < N; i++) {
            OvalRoi ovRoi = new OvalRoi(p[i][0]-p[i][2]+.5, p[i][1]-p[i][2]+.5, 2*p[i][2], 2*p[i][2]);
            ovRoi.setStrokeColor(c[lab[i]]);
            ovRoi.setName(Integer.toString(lab[i]));
            ov.add(ovRoi);
        }

        im.setOverlay(ov);
        im.show();

        ImageCanvas cnv = im.getCanvas();
        cnv.zoomIn(0, 0);
        cnv.zoomIn(0, 0);
        cnv.zoomIn(0, 0);
        cnv.zoomIn(0, 0);
        cnv.zoomIn(0, 0);

        // save it
        FileSaver fs = new FileSaver(im);
        fs.saveAsTiff("clustering2d.tif");

        /*
        example with indexes & inter-distances
         */

        IJ.log("#########################");
        IJ.log("#### TEST 2 ####");
        IJ.log("#########################");

        // random numbers in some interval, each has a value and index associated
        int N = 15;
        int RANGE = 40;
        int CLUSTER_TH = 1;
        int[] values = new int[N];
        int[] values_idx = new int[N];
        for (int i = 0; i < N; i++) {
            values[i] = new Random().nextInt(RANGE); // 0%(RANGE-1)
            values_idx[i] = i;
        }
        IJ.log("\nINPUT:\n");
        IJ.log(Arrays.toString(values));
        IJ.log("\nINDEXES:\n");
        IJ.log(Arrays.toString(values_idx));

        // table with precomputed inter-distances will be necessary for cluster mechanism
        int[][] dists = new int[N][N];
        for (int i=0; i< N; i++) {
            for (int j = i; j < N; j++) {
                if(i==j) {
                    dists[i][j] = 0;
                }
                else {
                    dists[i][j] = Math.abs(values[i]-values[j]);
                    dists[j][i] = Math.abs(values[i]-values[j]);
                }
            }
        }

        IJ.log("\nDISTS:\n");
        for (int k = 0; k < dists.length; k++) {
            IJ.log(Arrays.toString(dists[k]));
        }

        // show clusters in plot
        float[] plot_x = new float[values.length];
        for (int i=0; i<values.length; i++) plot_x[i] = (float) values[i];
        float[] plot_y = new float[values.length];

        Plot pl1 = new Plot("", "samples", "", new float[1], new float[1]);
        pl1.setLimits(0, RANGE, 0, 1.2);

        Arrays.fill(plot_y, 1);
        pl1.addPoints(plot_x, plot_y, Plot.CIRCLE);

        // plot vertical Lines for each value
        plot_y = new float[]{0, 1};
        plot_x = new float[2]; // {values[0], values[0]};
        for (int k=0; k<values.length; k++) {
            plot_x[0] = values[k];
            plot_x[1] = values[k];
            pl1.addPoints(plot_x, plot_y, Plot.LINE);
        }
        pl1.show();

        // cluster indexes of some values
        int[] lab2 = clustering(values_idx, dists, CLUSTER_TH); // lab2 will contain a cluster label for every value submitted
        // final outputs  - use labels to provide the final output
        ArrayList<float[]> ex = extracting(lab2, values); // there is a direct correspondence values_idx ~

        for (int i = 0; i < ex.size(); i++) {
            IJ.log("" + ex.get(i)[0] + " , " + ex.get(i)[1] + " elements");
        }

        Plot pl2 = new Plot("", "samples", "", new float[1], new float[1]);
        pl2.setLimits(0, RANGE, 0, 1.2);

        // plot vertical Lines for each value in it's cluster color

        c = new Color[N];
        for (int i = 0; i < N; i++) c[i] = getRandomColor();

        plot_y = new float[]{0, 1};
        plot_x = new float[2];//{values[0], values[0]};
        for (int k=0; k<values.length; k++) {
            plot_x[0] = values[k];
            plot_x[1] = values[k];
            pl2.setColor(c[lab2[k]]);
            pl2.setLineWidth(2);
            pl2.addPoints(plot_x, plot_y, Plot.LINE);
        }

        for (int k=0; k<values.length; k++) {
            pl2.setColor(c[lab2[k]]);
            pl2.addLabel((float) values[k]/RANGE, 0.1, Integer.toString(lab2[k]));

        }

        // show final clusters
        float[] cls_x = new float[ex.size()];
        float[] cls_y = new float[ex.size()];
        Arrays.fill(cls_y, 0.5f);
        for (int i=0; i<ex.size(); i++)
            cls_x[i] = ex.get(i)[0]; // /RANGE
        pl2.setColor(Color.BLACK);

        pl2.addPoints(cls_x, cls_y, Plot.TRIANGLE);

        pl2.show();

        IJ.log("done.");

    }
}
