package com.braincadet.phd.swc;

import ij.IJ;

import java.io.*;
import java.util.*;

/**
 * reader for swc reconstruction file
 *
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 *
 * revised for usage as a simple reader on 10/13/2015
 */
public class ReadSWC {

    // IDs have to be unique in swc format, therefore, IDs cannot repeat
    // depending on the export tool - it can happen that the IDs repeat or missing
    // input is read as the list of nodes (Node class) taking into account the linking
    // Node has the geometry + the link towards the neighbour

//    Standardized swc files (www.neuromorpho.org) -
//    0 - undefined
//    1 - soma
//    2 - axon
//    3 - (basal) dendrite
//    4 - apical dendrite
//    5+ - custom

    // new list with nodes
    public ArrayList<Node> nnodes = new ArrayList<Node>();

    public int maxID = Integer.MIN_VALUE;

    public float minR = Float.POSITIVE_INFINITY, maxR = Float.NEGATIVE_INFINITY;
    public float minX = Float.POSITIVE_INFINITY, maxX = Float.NEGATIVE_INFINITY;
    public float minY = Float.POSITIVE_INFINITY, maxY = Float.NEGATIVE_INFINITY;
    public float minZ = Float.POSITIVE_INFINITY, maxZ = Float.NEGATIVE_INFINITY;

    // indexes
    public static int ID 		= 0;
    public static int TYPE 		= 1;
    public static int XCOORD 	= 2;
    public static int YCOORD 	= 3;
    public static int ZCOORD 	= 4;
    public static int RADIUS 	= 5;
    public static int MOTHER 	= 6;

    public static int SWC_LINE_LENGTH = 7;

    public ReadSWC(String _swcFilePath) {

        String swcFilePath = new File(_swcFilePath).getAbsolutePath();// path to swc file

        if (!(new File(swcFilePath).exists())) {
            IJ.log(swcFilePath + " does not exist!");
            return;
        }

        IJ.log("reading...\t"+_swcFilePath);

        ArrayList<float[]> nodes_load = new ArrayList<float[]>(); // 1x7 rows (swc format)

        // read the node list of line elements - it's not guaranteed that the node IDs will be compatible with indexing
        // in sense that they are arranged ascending and even that all will exist alltogether
        // also can happen (although illegal to have doubling of the IDs)
        // read it first all to see what th efull range if IDs would be to know how much to allocate
        try { // scan the file

            FileInputStream fstream 	= new FileInputStream(swcFilePath);
            BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
            String read_line;

            while ( (read_line = br.readLine()) != null ) { // it will break on the empty line !!!

//              System.out.println("happened that it was empty ["+read_line+"]");//+br.readLine());//+"----->"+br.readLine()+"----"+(read_line = br.readLine()) != null);

                if (read_line.isEmpty()) continue;

                if(!read_line.trim().startsWith("#")) { // # are comments

//                    fileLength++;
                    // split values

                    String[] 	readLn = 	read_line.trim().replaceAll("," , ".").split("\\s+");

                    if (readLn.length!=SWC_LINE_LENGTH) continue; // skip the line that did not have enough values

                    float[] 	valsLn = 	new float[SWC_LINE_LENGTH]; // x, y, z, mother_index

                    valsLn[0] = Integer.valueOf(readLn[ID].trim()).intValue();      // id
                    valsLn[1] = Integer.valueOf(readLn[TYPE].trim()).intValue();    // type

                    valsLn[2] = Float.valueOf(readLn[XCOORD].trim()).floatValue();  // x, y, z
                    valsLn[3] = Float.valueOf(readLn[YCOORD].trim()).floatValue();
                    valsLn[4] = Float.valueOf(readLn[ZCOORD].trim()).floatValue();

                    valsLn[5] = Float.valueOf(readLn[RADIUS].trim()).floatValue();  // radius

                    valsLn[6] = Integer.valueOf(readLn[MOTHER].trim()).intValue();  // mother idx

                    nodes_load.add(valsLn);

                    maxID = (valsLn[0]>maxID)? (int) valsLn[0] : maxID;

                    minR = (valsLn[RADIUS]<minR)? valsLn[RADIUS] : minR;
                    maxR = (valsLn[RADIUS]>maxR)? valsLn[RADIUS] : maxR;

                    minX = (valsLn[XCOORD]<minX)? valsLn[XCOORD] : minX;
                    maxX = (valsLn[XCOORD]>maxX)? valsLn[XCOORD] : maxX;

                    minY = (valsLn[YCOORD]<minY)? valsLn[YCOORD] : minY;
                    maxY = (valsLn[YCOORD]>maxY)? valsLn[YCOORD] : maxY;

                    minZ = (valsLn[ZCOORD]<minZ)? valsLn[ZCOORD] : minZ;
                    maxZ = (valsLn[ZCOORD]>maxZ)? valsLn[ZCOORD] : maxZ;

                }

            }

            br.close();
            fstream.close();

        }
        catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }

        // analyze
        System.out.println(nodes_load.size() + " lines");
        System.out.println("maxID= " + maxID);

        // initialize all with null, ID will correspond to the index in nnodes list
        nnodes = new ArrayList<Node>(maxID + 1); // zeroth is dummy one (so that indexing can work with IDs) hence there is one more allocated
        for (int i = 0; i <= maxID; i++) nnodes.add(i, null);

        for (int i = 0; i < nodes_load.size(); i++) { // fill the nnodes list elements only

            int     currId      = Math.round(nodes_load.get(i)[ID]);
            int     currType    = Math.round(nodes_load.get(i)[TYPE]);

            float   currX       = nodes_load.get(i)[XCOORD];
            float   currY       = nodes_load.get(i)[YCOORD];
            float   currZ       = nodes_load.get(i)[ZCOORD];

            float   currR       = nodes_load.get(i)[RADIUS];

            int     prevId      = Math.round(nodes_load.get(i)[MOTHER]);

            if (nnodes.get(currId)==null) { // add the node
                nnodes.set(currId, new Node(currX, currY, currZ, currR, currType));
//                nnodes.get(i).nbr.add(prevId);
            }
            else {
                System.out.println("DOUBLING NODE ID HAPPENED!");
                // something is there, this means that we're doubling the same node (illegal but can happen)
                // by convention we'll assume that the doubled node is the same x,y,z,r - so we won't read the new coordinates again but the neighbours only
//                nnodes.get(i).nbr.add(prevId);
            }

        }

        // once the nodes are added add 2-directional connections read from the swc file (both neighbour references can be added as the nodes exist)
        for (int i = 0; i < nodes_load.size(); i++) {

            int     currId      = Math.round(nodes_load.get(i)[ID]);
            int     prevId      = Math.round(nodes_load.get(i)[MOTHER]);

            if (prevId!=-1) {
                nnodes.get(currId).nbr.add(prevId);
                nnodes.get(prevId).nbr.add(currId);
            }

        }

        // remove duplicate neighbouring links so that it does not confuse the bfs when checking connectivity
        remove_duplicate_neighbourhoods(nnodes);

    }

    private void remove_duplicate_neighbourhoods(ArrayList<Node> _nnodes) {

        // remove duplicate neighbourhood links
        System.out.print("\nremove duplicated links... ");
        for (int i = 0; i < _nnodes.size(); i++) {
            if (_nnodes.get(i)!=null) {
                Set<Integer> set = new HashSet<Integer>();
                set.addAll(_nnodes.get(i).nbr);
                _nnodes.get(i).nbr.clear();
                _nnodes.get(i).nbr.addAll(set);
            }
        }
        System.out.println("done.");

        // check if the neighborhoods are conistent
        System.out.print("\nchecking neighbourhood consistency... ");
        for (int i = 0; i < _nnodes.size(); i++) {
            if (_nnodes.get(i)!=null) {
                for (int j = 0; j < _nnodes.get(i).nbr.size(); j++) {
                    int nbr_idx = _nnodes.get(i).nbr.get(j);
                    if (nbr_idx>0) {
                        if (Collections.frequency(_nnodes.get(nbr_idx).nbr, i)!=1)
                            System.out.println("ERROR: " + i + " -- " + nbr_idx);
                    }
                    else {
//                        System.out.println("ERROR: parent id was " + nbr_idx);
                    }
                }
            }
        }
        System.out.println("done.");
    }

}