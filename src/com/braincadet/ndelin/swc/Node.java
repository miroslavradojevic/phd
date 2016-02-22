package com.braincadet.ndelin.swc;

import java.util.ArrayList;

/**
 * Created by miroslav on 13-3-15.
 */
public class Node {

    public float[]  loc;
    public float    r;
    public int      type;

    // types as in neuromorpho.org description
    public static int NOTHING = 0;
    public static int SOMA = 1;
    public static int AXON = 2;
    public static int BASAL_DENDRITE = 3;
    public static int APICAL_DENDRITE = 4;
    public static int CUSTOM1 = 5;
    public static int CUSTOM2 = 6;
    public static int UNDEFINED = 7;

    public ArrayList<Integer> nbr;

    public Node(float xn, float yn, float zn, float rn) {
        loc = new float[3];
        loc[0]  = xn;
        loc[1]  = yn;
        loc[2]  = zn;
        r       = rn;
        type    = NOTHING;
        nbr = new ArrayList<Integer>();
    }

    public Node(float xn, float yn, float zn, float rn, int typ) {
        this(xn, yn, zn, rn);
        type    = typ;
        nbr = new ArrayList<Integer>();
    }

//    public static boolean areNeighbours(int nidx1, ArrayList<Integer> nbrs1, int nidx2, ArrayList<Integer> nbrs2) {
//        return nbrs1.contains(nidx2) && nbrs2.contains(nidx1);
//    }

//    public static ArrayList<Integer> intersection(ArrayList<Integer> list1, ArrayList<Integer> list2) {
//
//        ArrayList<Integer> list = new ArrayList<Integer>();
//
//        for (Integer t : list1) {
//            if(list2.contains(t)) {
//                list.add(t);
//            }
//        }
//
//        return list;
//    }

//    public float overlap(Node node_to_check) { // [0,1]
//        // use volumetric overlap of two node spheres, normalized with the volume of the smaller sphere node
//        // http://mathworld.wolfram.com/Sphere-SphereIntersection.html
//
//        float d = (float) ((loc.length==2)? Math.sqrt(Math.pow(loc[0]-node_to_check.loc[0],2)+Math.pow(loc[1]-node_to_check.loc[1],2)) :
//                Math.sqrt(Math.pow(loc[0]-node_to_check.loc[0],2)+Math.pow(loc[1]-node_to_check.loc[1],2)+Math.pow(loc[2]-node_to_check.loc[2],2)));
//
//        float R1 = r;
//        float R2 = node_to_check.r;
//
//        if (d>R1+R2) return -1;
//
//        float V = (float) ((Math.PI/(12*d)) * Math.pow(R1 + R2 - d, 2) * (d*d+6*R1*R2    +2*d*R2-3*R2*R2     +2*d*R1-3*R1*R1)); // intersection volume using d, R1, R2
//
//        float Vnorm = (float) ((4/3f)*Math.PI*Math.pow(Math.min(R1,R2),3)); // (4/3)*pi*r^3 sphere volume
//
//        return (Vnorm>Float.MIN_VALUE)? (V/Vnorm) : 0f ; // return normalized volumetric overlap
//
//    }

}