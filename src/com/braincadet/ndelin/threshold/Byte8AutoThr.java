package com.braincadet.ndelin.threshold;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.process.*;

public class Byte8AutoThr implements PlugIn {

    private static final AutoThresholder.Method methods[] = {
            AutoThresholder.Method.Default,
//            AutoThresholder.Method.Huang,
//            AutoThresholder.Method.Intermodes,
            AutoThresholder.Method.IsoData,
//            AutoThresholder.Method.Li,
//            AutoThresholder.Method.MaxEntropy,
//            AutoThresholder.Method.Mean,
//            AutoThresholder.Method.MinError,
//            AutoThresholder.Method.Minimum,
            AutoThresholder.Method.Moments,
            AutoThresholder.Method.Otsu,
//            AutoThresholder.Method.Percentile,
//            AutoThresholder.Method.RenyiEntropy,
//            AutoThresholder.Method.Shanbhag,
            AutoThresholder.Method.Triangle,
//            AutoThresholder.Method.Yen
    };

    public void run(String s) {

        ImagePlus imp1 = IJ.getImage();

        if (!(imp1.getProcessor() instanceof ByteProcessor)) {
            IJ.log("Image needs to be BYTE8.");
            return;
        }

        int[] h = gethist(imp1, imp1.getStack().getSize()>1);

        if (h!=null) {

            AutoThresholder at = new AutoThresholder();

            for (int i = 0; i < methods.length; i++) {
                int threshold = at.getThreshold(methods[i], h);
                ImagePlus imp2 = imp1.duplicate();
                applythreshold(threshold, imp2); // will modify values
                imp2.setTitle(methods[i].toString());
                imp2.show();
                IJ.run(imp2, "Z Project...", "projection=[Max Intensity]");
                imp2.close();
            }

        }

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

    public int[] gethist(ImagePlus imp, boolean usestack) {

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
}