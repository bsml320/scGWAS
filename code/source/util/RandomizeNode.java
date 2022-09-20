package util;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.Random;
import java.util.Vector;

public class RandomizeNode {
	public static Hashtable<String, Double> randomize_node (Hashtable<String, Double> nodeHash){
		Object[] nodes_in_G    = nodeHash.keySet().toArray();
		Random r = new Random();
		double[] random_index_for_Y = new double[nodeHash.size()];
		double[] random_node_score = new double[nodeHash.size()];
		Hashtable<String, Double> random_nodeHash = new Hashtable<String, Double>();
		for(int i=0; i<nodeHash.size(); i++){
			random_index_for_Y[i] = r.nextDouble();
			random_node_score[i] = (double)nodeHash.get(nodes_in_G[i]);
		}
		
		// bubble sort random_index_for_Y, so that  random_Y is randomized
		boolean swapped = true; int jj = 0; double tmp = 0;
	    while (swapped) {
	        swapped = false;
	        jj++;
	        for (int i = 0; i < nodeHash.size() - jj; i++) {
	            if (random_index_for_Y[i] < random_index_for_Y[i + 1]) {
	                tmp = random_index_for_Y[i];
	                random_index_for_Y[i] = random_index_for_Y[i + 1];
	                random_index_for_Y[i + 1] = tmp;
	                
	                tmp = random_node_score[i];
	                random_node_score[i] = random_node_score[i + 1];
	                random_node_score[i + 1] = tmp;
	                
	                swapped = true;
	            }
	        }
	    } // end of bubble sort
	    
	    for(int k=0; k<nodes_in_G.length; k++){
	    	random_nodeHash.put((String)nodes_in_G[k], random_node_score[k]);
	    }
	    
		return(random_nodeHash);
	}
	
	public static Hashtable<String, Double> randomize_node_by_group (Hashtable<String, Double> nodeHash, Hashtable<String, Integer> groupHash){
		Hashtable<String, Double> random_nodeHash = new Hashtable<String, Double>();
		
		Object[] nodes_in_G    = groupHash.keySet().toArray();
		Hashtable<Integer, Vector<String>> group_label = new Hashtable<Integer, Vector<String>>();
		for(int i=0; i<groupHash.size(); i++){
			if(!nodeHash.containsKey(nodes_in_G[i]))continue;
			int group = groupHash.get(nodes_in_G[i]);
			if(group_label.containsKey(group)) {
				Vector<String> vec = group_label.get(group);
				vec.add( (String)nodes_in_G[i]);
				group_label.put(group, vec);
			} else {
				Vector<String> vec = new Vector<String>();
				vec.add( (String)nodes_in_G[i]);
				group_label.put(group, vec);
			}
		}
		
		Object[] unique_group_label = group_label.keySet().toArray();
		Arrays.sort(unique_group_label);
		for(int i=0; i<unique_group_label.length; i++) {
			int cur_group_label = (int) unique_group_label[i];
			Vector<String> vec = group_label.get(cur_group_label);
			int cur_group_size = vec.size();
			//System.out.println(cur_group_label+"\t"+cur_group_size);
			
			Random r = new Random();
			double[] random_index_for_Y = new double[cur_group_size];
			double[] random_node_score = new double[cur_group_size];
			for(int j=0; j<cur_group_size; j++){
				random_index_for_Y[j] = r.nextDouble();
				random_node_score[j] = (double)nodeHash.get(vec.get(j));
			}
			
			// bubble sort random_index_for_Y, so that  random_Y is randomized
			boolean swapped = true; int jj = 0; double tmp = 0;
		    while (swapped) {
		        swapped = false;
		        jj++;
		        for (int i1 = 0; i1 < cur_group_size - jj; i1++) {
		            if (random_index_for_Y[i1] < random_index_for_Y[i1 + 1]) {
		                tmp = random_index_for_Y[i1];
		                random_index_for_Y[i1] = random_index_for_Y[i1 + 1];
		                random_index_for_Y[i1 + 1] = tmp;
		                
		                tmp = random_node_score[i1];
		                random_node_score[i1] = random_node_score[i1 + 1];
		                random_node_score[i1 + 1] = tmp;
		                
		                swapped = true;
		            }
		        }
		    } // end of bubble sort
		    
		    for(int j=0; j<cur_group_size; j++){
		    	random_nodeHash.put((String)vec.get(j), random_node_score[j]);
		    }
		}
		
		return(random_nodeHash);
	}
	
	public static double[] randomize_edge (double[] edge_weights){
		int n_edges = edge_weights.length;
		Random r = new Random();
		double[] random_index_for_Y = new double[n_edges];
		double[] random_edge_score  = new double[n_edges];
		for(int i=0; i<n_edges; i++){
			random_index_for_Y[i] = r.nextDouble();
			random_edge_score[i] = edge_weights[i];
		}
		
		// bubble sort random_index_for_Y, so that  random_Y is randomized
		boolean swapped = true; int jj = 0; double tmp = 0;
	    while (swapped) {
	        swapped = false;
	        jj++;
	        for (int i = 0; i < n_edges - jj; i++) {
	            if (random_index_for_Y[i] < random_index_for_Y[i + 1]) {
	                tmp = random_index_for_Y[i];
	                random_index_for_Y[i] = random_index_for_Y[i + 1];
	                random_index_for_Y[i + 1] = tmp;
	                
	                tmp = random_edge_score[i];
	                random_edge_score[i] = random_edge_score[i + 1];
	                random_edge_score[i + 1] = tmp;
	                
	                swapped = true;
	            }
	        }
	    } // end of bubble sort
	    
		return(random_edge_score);
	}
	
}
