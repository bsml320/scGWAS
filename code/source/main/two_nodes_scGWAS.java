package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Random;
import java.util.Set;
import java.util.Vector;

import util.QuantileNormalization;
import util.RandomizeNode;
import util.rankify;

public class two_nodes_scGWAS {
	
	public static void scGWAS_wrapper(String gwas_node_weight_file, 
			                                   String scrn_weight_file, 
			                                   String module_score, 
			                                   String ppi_file, 
			                                   String normalization_model, 
			                                   Hashtable<String, String> cell_type_ht,
			                                   double r_include, 
			                                   double r_exclude, 
			                                   String outFile,
			                                   String int_genes_file,
			                                   String exclude_genes_file,
			                                   boolean permutation,
			                                   boolean remove_zero, 
			                                   int module_num,
			                                   boolean verbosity)throws Exception{
		
		String str;
		String[] line;
		
		// read GWAS node scores
		Hashtable<String, Double> nodeHash_gwas = new Hashtable<String, Double>();
		BufferedReader bin = new BufferedReader (new FileReader (new File(gwas_node_weight_file)));
		while((str=bin.readLine())!=null) {
			String[] tmp = str.split("\t");
			if(tmp[0].equalsIgnoreCase("Gene"))continue;
			nodeHash_gwas.put(tmp[0],Double.parseDouble(tmp[1]));
		}bin.close();
		System.out.println("Read in GWAS genes: "+nodeHash_gwas.size());
		///////////////////////////////////////////////////////////////////////////////////////////
		
		// average expression for each gene in each cell type
		Hashtable<String, Hashtable<String, Double>> sc_avgexpr_hash = new Hashtable<String, Hashtable<String, Double>>();
		Hashtable<String, Integer> sc_header_index = new Hashtable<String, Integer>();
		bin = new BufferedReader (new FileReader (new File(scrn_weight_file)));
		str=bin.readLine(); // skip the header line.
		line = str.split("\t");
		for(int i=1; i<line.length; i++) {
			if(cell_type_ht!=null) { // read only input cell types
				if(cell_type_ht.containsKey(line[i])) { // consider only input cell types that are also in the header line
					sc_header_index.put(line[i], i);
					sc_avgexpr_hash.put(line[i], new Hashtable<String, Double>());
				}
			} else { // read in all cell types
				sc_header_index.put(line[i], i);
				sc_avgexpr_hash.put(line[i], new Hashtable<String, Double>());
			}
		}
		System.out.println("Identified " + sc_avgexpr_hash.size() + " cell types for the following analyses");
		///////////////////////////////////////////////////////////////////////////////////////////
		
		Object[] cell_type_mat = sc_header_index.keySet().toArray();
		Arrays.sort(cell_type_mat);
		Hashtable<String, String> raw_scrn_genes = new Hashtable<String, String>();
		while((str=bin.readLine())!=null) {
			line = str.split("\t");
			raw_scrn_genes.put(line[0], "");
			
			for(int i=0; i<cell_type_mat.length; i++) {
				int index = (int)sc_header_index.get(cell_type_mat[i]);
				double avg_expr = Double.parseDouble(line[index]);
				Hashtable<String, Double> nodeHash_sc = (Hashtable<String, Double>)sc_avgexpr_hash.get(cell_type_mat[i]);
				if(remove_zero) {
					if(avg_expr < 2.2e-16)continue;
				}
				nodeHash_sc.put(line[0], avg_expr);
				sc_avgexpr_hash.put((String) cell_type_mat[i], nodeHash_sc);
			}
			
		}bin.close();
		
		///////////////////////////////////////////////////////////////////////////////////////////
		Hashtable<String, String> intGenes = new Hashtable<String, String>();
		if(int_genes_file!=null) {
			File int_genes_file_F = new File(int_genes_file);
			if(!int_genes_file_F.exists()) {
				System.out.println("File "+int_genes_file+" does not exist. Failed to read in the genes of interest. The following analsis would be focusing on all genes.");
			} else {
				bin = new BufferedReader (new FileReader (int_genes_file_F));
				while((str=bin.readLine())!=null) {
					intGenes.put(str, "");
				}
				bin.close();
			}
			System.out.println("Read in n = "+intGenes.size()+" genes of interest. The following analysis would be only focusing on these genes.");
		} else {
			intGenes = null;
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		
		///////////////////////////////////////////////////////////////////////////////////////////
		Hashtable<String, String> excludeGenes = new Hashtable<String, String>();
		if (exclude_genes_file != null) {
			File exclude_genes_file_F = new File(exclude_genes_file);
			if (!exclude_genes_file_F.exists()) {
				System.out.println("File " + exclude_genes_file	+ " does not exist.");
			} else {
				bin = new BufferedReader(new FileReader(exclude_genes_file_F));
				while ((str = bin.readLine()) != null) {
					excludeGenes.put(str, "");
				}
				bin.close();
			}
			System.out.println("Read in n = " + excludeGenes.size()	+ " genes to be excluded. ");
		} else {
			excludeGenes = null;
		}
		///////////////////////////////////////////////////////////////////////////////////////////

		Vector<String[]> PPI = new Vector<String[]>();
		Hashtable<String, Vector<Integer>> gene2PPI = new Hashtable<String, Vector<Integer>>();
		int idx = 0;
		bin = new BufferedReader (new FileReader (new File(ppi_file)));
		int count_PPI = 0;
		while((str=bin.readLine())!=null){
			String[] tmp = str.split("\t");
			count_PPI++;
			
			if(!nodeHash_gwas.containsKey(tmp[0]) | !nodeHash_gwas.containsKey(tmp[1]) )continue;
			if(!raw_scrn_genes.containsKey(tmp[0]) | !raw_scrn_genes.containsKey(tmp[1]) )continue;
			
			if(excludeGenes!=null) {
				if(excludeGenes.containsKey(tmp[0]) | excludeGenes.containsKey(tmp[1]) )continue;
			}
			
			PPI.add(tmp);
			if(gene2PPI.containsKey(tmp[0])){
				Vector<Integer> idxVec = (Vector<Integer>)gene2PPI.get(tmp[0]);
				idxVec.add(idx);
				gene2PPI.put(tmp[0], idxVec);
			} else {
				Vector<Integer> idxVec = new Vector<Integer>();
				idxVec.add(idx);
				gene2PPI.put(tmp[0], idxVec);
			}
			if(gene2PPI.containsKey(tmp[1])){
				Vector<Integer> idxVec = (Vector<Integer>)gene2PPI.get(tmp[1]);
				idxVec.add(idx);
				gene2PPI.put(tmp[1], idxVec);
			} else {
				Vector<Integer> idxVec = new Vector<Integer>();
				idxVec.add(idx);
				gene2PPI.put(tmp[1], idxVec);
			}
			idx = idx + 1;
		}bin.close(); 
		System.out.println("# input GWAS genes = "+nodeHash_gwas.size() + 
				"\n# scRNAseq genes = " + raw_scrn_genes.size() +"\n# raw input PPI = "+count_PPI+
				"\nWorking data set:\n# PPIs = "+PPI.size()+" edges (both nodes required to be present in the GWAS file and the scRNAseq file)\n# nodes = "+gene2PPI.size()+
				" (with weights from both GWAS and scRNAseq files)");
		
		///////////////////////////////////////////////////////////////////////////////////////////
		// define gene group by degree
		Object[] genes_in_G   = gene2PPI.keySet().toArray();
		int[] degree = new int[genes_in_G.length];
		Hashtable<String, Integer> groupHash  = new Hashtable<String, Integer>();
		for(int i=0; i<genes_in_G.length; i++) {
			degree[i] = gene2PPI.get(genes_in_G[i]).size();
		}
		
		boolean swapped = true;
	    int bubble_jj = 0;
	    int tmp_degree; Object tmp_gene;
	    while (swapped) {
	        swapped = false;
	        bubble_jj++;
	        for (int i = 0; i < degree.length - bubble_jj; i++) {
	            if (degree[i] < degree[i + 1]) {
	            	tmp_degree = degree[i];
	            	degree[i] = degree[i + 1];
	            	degree[i + 1] = tmp_degree;
	            	
	            	tmp_gene = genes_in_G[i];
	            	genes_in_G[i] = genes_in_G[i + 1];
	            	genes_in_G[i + 1] = tmp_gene;
	                swapped = true;
	            }
	        }
	    }
	    
	    int group_size = (int)Math.ceil(genes_in_G.length/100.0);
		System.out.println("Bin size: "+group_size);
		for(int i=0; i<100; i++) {
			int start_idx = i * group_size, end_idx = (i+1) * group_size;
			if(end_idx >= genes_in_G.length)end_idx = genes_in_G.length;
			for(int j=start_idx; j<end_idx; j++) {
				groupHash.put( (String)genes_in_G[j], i);
			}
		}
		
		///////////////////////////////////////////////////////////////////////////////////////////
		
		Hashtable<String, Hashtable<String, Double>> header_new_nodeHash_gwas = new Hashtable<String, Hashtable<String, Double>>();
		Hashtable<String, Hashtable<String, Double>> header_new_nodeHash_scrn = new Hashtable<String, Hashtable<String, Double>>();
		for(int i=0; i<cell_type_mat.length; i++) {
			String header_name = (String)cell_type_mat[i];
			Hashtable<String, Double> nodeHash_sc = (Hashtable<String, Double>)sc_avgexpr_hash.get(header_name);
			
			File f = new File(outFile+"_"+header_name+".modules.txt");
			if(f.exists()) {
				System.out.println("Cell type: " + header_name+" already exist.");
				continue;
			}
			
			Date dNow = new Date();
			SimpleDateFormat ft = new SimpleDateFormat("E yyyy.MM.dd 'at' hh:mm:ss a zzz");
			System.out.println("\nStart working on column <"+header_name+"> at "+ft.format(dNow));
			System.out.println("Current scRNAseq data: genes = "+nodeHash_sc.size());
			
			// normalization: calibration
			Vector<Hashtable<String, Double>> vec = new Vector<Hashtable<String, Double>>();
			Hashtable<String, Double> new_nodeHash_gwas = new Hashtable<String, Double>();
			Hashtable<String, Double> new_nodeHash_sc = new Hashtable<String, Double>();
			
			System.out.println("normalization_model = "+normalization_model);
			if(normalization_model == "calibration") {
				vec = get_normalization_calibration(nodeHash_gwas, nodeHash_sc, verbosity, outFile, header_name, i);
				new_nodeHash_gwas = (Hashtable<String, Double>)vec.get(0);
				new_nodeHash_sc = (Hashtable<String, Double>)vec.get(1);
			} else {
				new_nodeHash_gwas = nodeHash_gwas;
				new_nodeHash_sc = nodeHash_sc;
			}
			
			System.out.println("new_nodeHash_gwas = "+new_nodeHash_gwas.size());
			System.out.println("new_nodeHash_sc   = "+new_nodeHash_sc.size());
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// THE scGWAS
			Hashtable<String, String> modules = scGWAS_two_nodes(new_nodeHash_gwas, new_nodeHash_sc, PPI, gene2PPI, r_include, r_exclude, module_score, intGenes);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			header_new_nodeHash_gwas.put(header_name, new_nodeHash_gwas);
			header_new_nodeHash_scrn.put(header_name, new_nodeHash_sc);
			
			System.out.println("# find modules = "+modules.size());
			BufferedWriter bout = new BufferedWriter(new FileWriter(new File(outFile+"_"+header_name+".modules.txt")));
			Object[] genes_in   = modules.keySet().toArray();
			Hashtable<String, String> module_size   = new Hashtable<String, String>();
			int max_module_size = 0;
			
			bout.write("module_genes\tseed\tedges\tmodule_score\tz_gwas\tz_scrna\n");
			for(int k=0; k<genes_in.length; k++){
				String tmp = (String)modules.get((String)genes_in[k]);
				String[] nodes = ((String)genes_in[k]).split(":");
				if(nodes.length > max_module_size) max_module_size = nodes.length;
				
				if(module_size.containsKey(nodes.length+"")) {
					String str_size = (String)module_size.get(nodes.length+"");
					double n_modules_at_size = Double.parseDouble(str_size) + 1;
					module_size.put(nodes.length+"", n_modules_at_size+"");
				} else {
					module_size.put(nodes.length+"", "1");
				}
				
				bout.write(genes_in[k]+"\t"+tmp+"\n");
			}
			bout.flush(); bout.close();
			
			// n_genes: the size for qualified modules to be tested by permutation
			// all_sizes: all possible sizes observed in real; extreme sizes such as 2 or those with <20 modules would not be tested by permutation
			
			Hashtable<String, double[]> n_genes   = new Hashtable<String, double[]>();
			Object[] all_sizes   = module_size.keySet().toArray();
			for(int k=0; k<all_sizes.length; k++) {
				String str_size = (String)all_sizes[k];
				double n_modules_at_size = Double.parseDouble(module_size.get(str_size));
				System.out.println("True modules:\t"+str_size+"\t"+n_modules_at_size);
				
				
				if(Double.parseDouble(str_size) == 2)continue;
				//if(n_modules_at_size < 20)continue; // in future, consider use this.
				
				if(Double.parseDouble(str_size) > 3 & n_modules_at_size < 20)continue;
				n_genes.put(str_size, new double[6]);
			}
			
			if(permutation){
				Object[] n_genes_per_module   = n_genes.keySet().toArray();
				Arrays.sort(n_genes_per_module);
				
				// perform virtual dmGWAS but using random node scores
				bout = new BufferedWriter(new FileWriter(new File(outFile+"_"+header_name+".random_modules.txt")));
				Hashtable<String, Vector<Double[]>> random_modules_by_size   = new Hashtable<String, Vector<Double[]>>();
				int round = 0;
				while(round < 100){
					Hashtable<String, Double> random_nodeHash_gwas = null;
					Hashtable<String, Double> random_nodeHash_scrn = null;
					
					random_nodeHash_gwas = RandomizeNode.randomize_node(new_nodeHash_gwas);
					random_nodeHash_scrn = RandomizeNode.randomize_node_by_group(new_nodeHash_sc, groupHash);
					
					Hashtable<String, String> random_modules = scGWAS_two_nodes(random_nodeHash_gwas, random_nodeHash_scrn, PPI, gene2PPI, r_include, r_exclude, module_score, null);
				    
				    Object[] random_genes_in = random_modules.keySet().toArray();
					for(int k=0; k<random_genes_in.length; k++){
						String tmp = (String)random_modules.get((String)random_genes_in[k]);
						String[] nodes = ((String)random_genes_in[k]).split(":");
						String[] module_info = tmp.split("\t");
						String size = nodes.length+"";
						if(random_modules_by_size.containsKey(size)) {
							Vector<Double[]> already_modules = random_modules_by_size.get(size);
							double m = Double.parseDouble(module_info[2]);
							double m_gwas = Double.parseDouble(module_info[3]);
							double m_scrn = Double.parseDouble(module_info[4]);
							already_modules.add(new Double[] {m, m_gwas, m_scrn});
							random_modules_by_size.put(size, already_modules);
						} else {
							Vector<Double[]> already_modules = new Vector<Double[]>();
							double m = Double.parseDouble(module_info[2]);
							double m_gwas = Double.parseDouble(module_info[3]);
							double m_scrn = Double.parseDouble(module_info[4]);
							already_modules.add(new Double[] {m, m_gwas, m_scrn});
							random_modules_by_size.put(size, already_modules);
						}
						bout.write(random_genes_in[k]+"\t"+tmp+"\t"+round+"\n");
					}
					bout.flush();
					round++;
					System.out.println();
					
					int boo_size_check = 0;
					for(int rr=0; rr<n_genes_per_module.length; rr++) {
						String size_str = (String)n_genes_per_module[rr];
						
						if(random_modules_by_size.containsKey(size_str)) {
							Vector<Double[]> already_modules = random_modules_by_size.get(size_str);
							
							if( Integer.parseInt(size_str) == max_module_size || size_str.equals("2") ) {
								if(already_modules.size() < module_num & round < 50) {
									boo_size_check++;
								}
								System.out.println("round = "+round+"\t"+size_str+"\t"+already_modules.size()+"\t"+boo_size_check+"\t(max size)");
							} else {
								if(already_modules.size() < module_num) {
									boo_size_check++;
								}
								System.out.println("round = "+round+"\t"+size_str+"\t"+already_modules.size()+"\t"+boo_size_check);
							}
							
						} else {
							boo_size_check++;
							System.out.println(round+"\t"+size_str+"\tmissing\t"+boo_size_check);
						} 
					}
					
					if(boo_size_check==0)break;
					
				}bout.flush(); bout.close();
				
				// get the mean and sd for each module size
				for(int rr=0; rr<n_genes_per_module.length; rr++) {
					String size_str = (String)n_genes_per_module[rr];
					if(!random_modules_by_size.containsKey(size_str)) {
						double size_info[] = new double[] {0, 0, 0, 0, 0, 0};
						n_genes.put(size_str, size_info);
						continue;
					}
					
					Vector<Double[]> already_modules = random_modules_by_size.get(size_str);
					double m=0, sd=0, gwas_mean = 0, gwas_sd = 0, scrn_mean = 0, scrn_sd = 0;
					for(int k1=0; k1<already_modules.size(); k1++){
						Double[] m_scores = already_modules.get(k1);
						m = m + m_scores[0];
						gwas_mean = gwas_mean + m_scores[1];
						scrn_mean = scrn_mean + m_scores[2];
					}
					m = m/already_modules.size();
					gwas_mean = gwas_mean/already_modules.size();
					scrn_mean = scrn_mean/already_modules.size();
					
					for(int k1=0; k1<already_modules.size(); k1++){
						Double[] m_scores = already_modules.get(k1);
					    sd = sd + Math.pow(m_scores[0] - m, 2);
					    gwas_sd = gwas_sd + Math.pow(m_scores[1] - gwas_mean, 2);
					    scrn_sd = scrn_sd + Math.pow(m_scores[2] - scrn_mean, 2);
					}
					sd = Math.sqrt(sd/(already_modules.size()-1) );
					gwas_sd = Math.sqrt(gwas_sd/(already_modules.size()-1) );
					scrn_sd = Math.sqrt(scrn_sd/(already_modules.size()-1) );
					
					double size_info[] = new double[] {m, sd, gwas_mean, gwas_sd, scrn_mean, scrn_sd};
					if(already_modules.size() < 5) {
						m=0; sd=0; gwas_mean = 0; gwas_sd = 0; scrn_mean = 0; scrn_sd = 0;
						size_info = new double[] {0, 0, gwas_mean, gwas_sd, scrn_mean, scrn_sd};
					}
					n_genes.put(size_str, size_info);
					System.out.println(size_str+"\t"+m+"\t"+sd+"\t"+gwas_mean+"\t"+gwas_sd+"\t"+scrn_mean+"\t"+scrn_sd);
				}
				
				// calculate module score based on size distribution
				for(int k=0; k<genes_in.length; k++){
					String[] nodes = ((String)genes_in[k]).split(":");  // module genes
					String size = nodes.length + "";
					
					if(!n_genes.containsKey(size))continue;
					
					String module_nodes_str = (String)genes_in[k]; 
			    	String info = (String)modules.get(module_nodes_str); // module info
			    	String[] info_array = info.split("\t");
			    	
			    	
			    	double size_info[] = n_genes.get(size);
			    	if(size_info[0] == 0 & size_info[1] == 0)continue;
			    	
			    	double z_score = (Double.parseDouble(info_array[2]) - size_info[0] )/size_info[1];
					double z_gwas  = (Double.parseDouble(info_array[3]) - size_info[2] )/size_info[3];
					double z_scrn  = (Double.parseDouble(info_array[4]) - size_info[4] )/size_info[5];
					info = info+"\t"+z_score+"\t"+z_gwas+"\t"+z_scrn;
					modules.put(module_nodes_str, info);
				}
				
				bout = new BufferedWriter(new FileWriter(new File(outFile+"_"+header_name+".modules.txt")));
				bout.write("module_genes\tseed\tedges\tmodule_score\tm_gwas\tm_scrnaseq\tmodule_score_z\tz_gwas\tz_scrnaseq\n");
				for(int k=0; k<genes_in.length; k++){
			    	String module_nodes_str = (String)genes_in[k];
			    	
			    	//String[] nodes = ((String)genes_in[k]).split(":");
					//if(nodes.length>=7 | nodes.length == 2)continue;
					
			    	String info = (String)modules.get(module_nodes_str);
			    	bout.write(module_nodes_str+"\t"+info+"\n");
				}
				bout.flush(); bout.close();
			}
			
		}
		
		
		BufferedWriter bout = new BufferedWriter(new FileWriter(new File(outFile+".combined-z.txt")));
		if(permutation) {
			bout.write("module_genes\tseed\tedges\tmodule_score\tm_gwas\tm_scrnaseq\tmodule_score_z\tz_gwas\tz_scrnaseq\tcell_type\n");
		} else {
			bout.write("module_genes\tseed\tedges\tmodule_score\tm_gwas\tm_scrnaseq\tmodule_score_z\tz_gwas\n");
		}
		for(int i=0; i<cell_type_mat.length; i++) {
			String header_name = (String)cell_type_mat[i];
			File f1 = new File(outFile+"_"+header_name+".modules.txt");
			if(!f1.exists()) {
				System.out.println("File not collected: "+f1.getAbsolutePath()+"/"+f1.getName());
				continue;
			}
			
			bin = new BufferedReader(new FileReader(new File(outFile+"_"+header_name+".modules.txt")));
			str=bin.readLine();
			while((str=bin.readLine())!=null){
				bout.write(str+"\t"+header_name+"\n");
			}bin.close();
			
			f1.delete();
		}bout.flush(); bout.close(); 
		
		if(permutation) {
			bout = new BufferedWriter(new FileWriter(new File(outFile+".all.random_modules.txt")));
			for(int i=0; i<cell_type_mat.length; i++) {
				String header_name = (String)cell_type_mat[i];
				
				File f1 = new File(outFile+"_"+header_name+".random_modules.txt");
				if(!f1.exists()) {
					System.out.println("File not collected: "+f1.getAbsolutePath()+"/"+f1.getName());
					continue;
				}
				
				bin = new BufferedReader(new FileReader(f1));
				str=bin.readLine();
				while((str=bin.readLine())!=null){
					bout.write(str+"\t"+header_name+"\n");
				}bin.close();
				
				f1.delete();
			}bout.flush(); bout.close();
		}
		
	}
	
	
	public static Hashtable<String, String> scGWAS_two_nodes (Hashtable<String, Double> nodeHash_gwas, 
			                                                  Hashtable<String, Double> nodeHash_sc, 
			                                                  Vector<String[]> PPI, 
			                                                  Hashtable<String, Vector<Integer>> gene2PPI, 
			                                                  double r_include, 
			                                                  double r_exclude, 
			                                                  String module_score, 
			                                                  Hashtable<String, String> intGenes){
		
		Hashtable<String, String> modules = new Hashtable<String, String>();
		String max_gene = "";
		double count = 1;
		
		Object[] genes = nodeHash_gwas.keySet().toArray();
		for(int j=0; j<genes.length; j++){
			
			if( Integer.parseInt(""+10*j/genes.length) == count){
				System.out.print( (count*10 ) +"%...");
				count++;
			}
			
			String seed_gene = (String)genes[j];
			if(!nodeHash_sc.containsKey(seed_gene) | !gene2PPI.containsKey(seed_gene))continue;
			if(intGenes!=null) {
				if(!intGenes.containsKey(seed_gene))continue;
			}
			
			//if(!seed_gene.equals("TNF"))continue;
			
			Hashtable<String, String> M = new Hashtable<String, String>();
			
			double[] new_score_best = new double[3];
			M.put(seed_gene, seed_gene);
			while(true){
				Object[] genes_in = M.keySet().toArray();
				max_gene = "";
				
				/*
				System.out.print("Current genes: ");
				for(int n=0; n<genes_in.length; n++){
					String gene = (String)genes_in[n];
					System.out.print(gene+",");
				}System.out.print("\n");
				*/
				
				// obtain all interactors
				Set<Integer> set = new HashSet<Integer>(new Vector<Integer>());
				for(int n=0; n<genes_in.length; n++){
					String gene = (String)genes_in[n];
					
					Vector<Integer> cur_idxVec = (Vector<Integer>)gene2PPI.get(gene);
					set.addAll(cur_idxVec);
				}
				
				double[] cur_score_1 = get_module_score(nodeHash_gwas, nodeHash_sc, genes_in, module_score);
				//System.out.print("\tcurrent_score_m = "+cur_score_1[0]+"\tm_gwas = "+cur_score_1[1]+"\tm_scrna = "+cur_score_1[2]+"\n");
				
				// find the best node in neighbor
				Hashtable<String, String> can_nodes = new Hashtable<String, String>();
				for(Integer k:set){
					String[] tmp_edge = (String[]) PPI.get(k);
					if(nodeHash_sc.containsKey(tmp_edge[0]) && nodeHash_sc.containsKey(tmp_edge[1])) {
						can_nodes.put(tmp_edge[0], "");
						can_nodes.put(tmp_edge[1], "");
					}
				}
				
				double current_n_genes = genes_in.length;
				new_score_best = cur_score_1;
				Object[] can_genes = can_nodes.keySet().toArray();
				for(int n1=0; n1<can_genes.length; n1++){
					String gene = (String)can_genes[n1];
					if(M.containsKey(gene))continue;
					
					double tmp_score_gwas = ( cur_score_1[1] * Math.sqrt(current_n_genes) + nodeHash_gwas.get(gene) )/Math.sqrt(current_n_genes+1);
					double tmp_score_scrn = ( cur_score_1[2] * Math.sqrt(current_n_genes) + nodeHash_sc.get(gene)   )/Math.sqrt(current_n_genes+1);
					double m = (tmp_score_gwas + tmp_score_scrn)/2;
					double tmp_module_score = tmp_score_gwas + tmp_score_scrn - Math.sqrt( (tmp_score_gwas - m)*(tmp_score_gwas - m) + (tmp_score_scrn - m)*(tmp_score_scrn - m) );
					
					double[] tmp_score_1 = new double[] {tmp_module_score, tmp_score_gwas, tmp_score_scrn};
					if(tmp_score_1[0] > new_score_best[0] & tmp_score_1[1] > new_score_best[1] & tmp_score_1[2] > new_score_best[2]){
						new_score_best = tmp_score_1; 
						max_gene = gene;
					}
					
				}
				
				if(new_score_best[0] > cur_score_1[0] * (1+r_include) & !max_gene.equals("") ){
					M.put(max_gene, ""); 
				} else {	
					break;	
				}
				
				//////////////////////////////////////////////////////////////////////////
				// exclude edge
				while(M.size() >= 3) {
					boolean boo = false;
					genes_in = M.keySet().toArray();
					current_n_genes = genes_in.length;
					
					Vector<String> nodeVec = get_vector(genes_in);
					cur_score_1 = get_module_score(nodeHash_gwas, nodeHash_sc, nodeVec, module_score);
					
					Set<Integer> set_in_module = get_edge_set(gene2PPI, PPI, nodeVec); 
					
					// obtain the max score after removing each edge
					String to_remove_gene_1 = "", to_remove_gene_2 = "";
					Hashtable<String, Double> edge_record = new Hashtable<String, Double>();
					for (Integer k : set_in_module) {
						String[] edge = (String[]) PPI.get(k);
						
						double tmp_score_gwas = ( cur_score_1[1] * Math.sqrt(current_n_genes) - nodeHash_gwas.get(edge[0]) - nodeHash_gwas.get(edge[1]) )/Math.sqrt(current_n_genes-2);
						double tmp_score_scrn = ( cur_score_1[2] * Math.sqrt(current_n_genes) - nodeHash_sc.get(edge[0])   - nodeHash_sc.get(edge[1])   )/Math.sqrt(current_n_genes-2);
						double me = (tmp_score_gwas + tmp_score_scrn)/2;
						double tmp_module_score = tmp_score_gwas + tmp_score_scrn - Math.sqrt( (tmp_score_gwas - me)*(tmp_score_gwas - me) + (tmp_score_scrn - me)*(tmp_score_scrn - me) );
						double[] tmp_score_1 = new double[] {tmp_module_score, tmp_score_gwas, tmp_score_scrn};
						if(tmp_score_1[0] > cur_score_1[0] * (1-r_exclude) ){
									edge_record.put(edge[0]+":"+edge[1], tmp_score_1[0]);
						}
					}

					if (edge_record.size() > 0) {
						
						// order edges
						Object[] candidate_remove_genes = (Object[])edge_record.keySet().toArray();
						double[] candidate_remove_genes_score = new double[candidate_remove_genes.length];
						for(int i=0; i<candidate_remove_genes.length; i++){
							candidate_remove_genes_score[i] = (Double)edge_record.get((String)candidate_remove_genes[i]);
							
						}
						
						// bubble sort
						boolean swapped = true;
					    int jj = 0;
					    double tmp; String tmp2="";
					    while (swapped) {
					        swapped = false;
					        jj++;
					        for (int i = 0; i < candidate_remove_genes_score.length - jj; i++) {
					            if (candidate_remove_genes_score[i] < candidate_remove_genes_score[i + 1]) {
					                tmp = candidate_remove_genes_score[i];
					                candidate_remove_genes_score[i] = candidate_remove_genes_score[i + 1];
					                candidate_remove_genes_score[i + 1] = tmp;
					                tmp2 = (String)candidate_remove_genes[i];
					                candidate_remove_genes[i] = candidate_remove_genes[i + 1];
					                candidate_remove_genes[i + 1] = tmp2;
					                swapped = true;
					            }
					        }
					    }
					    // end of ordering edges
					    
						for(int n=0; n<edge_record.size(); n++){
							String to_remove_gene = (String)candidate_remove_genes[n];
							to_remove_gene_1 = to_remove_gene.split(":")[0];
							to_remove_gene_2 = to_remove_gene.split(":")[1];
							
							// check if it's still connected
							genes_in = M.keySet().toArray();
							Vector<String> new_genes_in_vec = new Vector<String>();
							for (int m2 = 0; m2 < genes_in.length; m2++) {
								if (genes_in[m2].equals(to_remove_gene_1) || genes_in[m2].equals(to_remove_gene_2)) continue;
								new_genes_in_vec.add((String)genes_in[m2]);
							}
							set = get_edge_set(gene2PPI, PPI, new_genes_in_vec); 
							
							int[][] remove_edge_Edges = new int[set.size()][2];
							int kk = 0;
							for (Integer k : set) {
								String[] tmp_edge = (String[]) PPI.get(k);
								int from = -1, to = -1;
								for (int m2 = 0; m2 < new_genes_in_vec.size(); m2++) {
									String gene = (String) new_genes_in_vec.get(m2);
									if (gene.equals(to_remove_gene_1) || gene.equals(to_remove_gene_2))
										continue;
									if (tmp_edge[0].equals(gene)) {
										from = m2;
									}
									if (tmp_edge[1].equals(gene)) {
										to = m2;
									}
								}
								if (from != -1 && to != -1) {
									remove_edge_Edges[kk][0] = from;
									remove_edge_Edges[kk][1] = to;
									kk = kk + 1;
								}
							}
							int N_com = countComponents(new_genes_in_vec.size(), remove_edge_Edges);
							if (N_com != 1) continue;  // not connected

							M.remove(to_remove_gene_1); 
							M.remove(to_remove_gene_2); 
							
							//System.out.println("Exclude: "+to_remove_gene_1 + ":" + to_remove_gene_2);
							
							boo = true;
							break;
						}
						if(!boo)break;
					} else {
						break;
					}
				}
				// End of edge remove
				////////////////////////////////////////////////////////////////////////////////////////
				
				
				
				////////////////////////////////////////////////////////////////////////////////////////
				// node remove
				while (M.size() > 1) {
					boolean boo = false;

					// current module
					genes_in = M.keySet().toArray();
					cur_score_1 = get_module_score(nodeHash_gwas, nodeHash_sc, genes_in, module_score);
					current_n_genes = genes_in.length;
					
					// check each gene, record the scores should the gene be removed
					Hashtable<String, Double> record = new Hashtable<String, Double>();
					for (int m = 0; m < genes_in.length; m++) {
						String gene = (String) genes_in[m];
						if (gene.equals(max_gene))
							continue;
						
						double tmp_score_gwas = ( cur_score_1[1] * Math.sqrt(current_n_genes) - nodeHash_gwas.get(gene) )/Math.sqrt(current_n_genes-1);
						double tmp_score_scrn = ( cur_score_1[2] * Math.sqrt(current_n_genes) - nodeHash_sc.get(gene)   )/Math.sqrt(current_n_genes-1);
						double me = (tmp_score_gwas + tmp_score_scrn)/2;
						double tmp_module_score = tmp_score_gwas + tmp_score_scrn - Math.sqrt( (tmp_score_gwas - me)*(tmp_score_gwas - me) + (tmp_score_scrn - me)*(tmp_score_scrn - me) );
						double[] tmp_score_1 = new double[] {tmp_module_score, tmp_score_gwas, tmp_score_scrn};
						if (tmp_score_1[0] > cur_score_1[0] * (1 - r_exclude)) {
							record.put(gene, tmp_score_1[0]);
						}
					}

					if (record.size() > 0) {
						Object[] candidate_remove_genes = (Object[]) record.keySet().toArray();
						double[] candidate_remove_genes_score = new double[candidate_remove_genes.length];
						for (int i = 0; i < candidate_remove_genes.length; i++) {
							candidate_remove_genes_score[i] = (Double) record.get((String) candidate_remove_genes[i]);
						}

						// bubble sort
						boolean swapped = true;
						int bubble_jj = 0;
						double tmp;
						String tmp2 = "";
						while (swapped) {
							swapped = false;
							bubble_jj++;
							for (int i = 0; i < candidate_remove_genes_score.length - bubble_jj; i++) {
								if (candidate_remove_genes_score[i] < candidate_remove_genes_score[i + 1]) {
									tmp = candidate_remove_genes_score[i];
									candidate_remove_genes_score[i] = candidate_remove_genes_score[i + 1];
									candidate_remove_genes_score[i + 1] = tmp;
									tmp2 = (String) candidate_remove_genes[i];
									candidate_remove_genes[i] = candidate_remove_genes[i + 1];
									candidate_remove_genes[i + 1] = tmp2;
									swapped = true;
								}
							}
						}

						for (int i = 0; i < candidate_remove_genes.length; i++) {
							String to_remove_gene = (String) candidate_remove_genes[i];

							/////////////////////////////////////////////////// check if it's still connected
							genes_in = M.keySet().toArray();
							Vector<String> new_genes_in_vec = new Vector<String>();
							for (int m2 = 0; m2 < genes_in.length; m2++) {
								if (genes_in[m2].equals(to_remove_gene))
									continue;
								new_genes_in_vec.add((String) genes_in[m2]);
							}

							set = get_edge_set(gene2PPI, PPI, new_genes_in_vec);
							
							int[][] remove_node_Edges = new int[set.size()][2];
							int kk = 0;
							for (Integer k : set) {
								String[] tmp_edge = (String[]) PPI.get(k);
								int from = -1, to = -1;
								for (int m2 = 0; m2 < new_genes_in_vec.size(); m2++) {
									String gene = (String) new_genes_in_vec.get(m2);
									if (gene.equals(to_remove_gene))
										continue;
									if (tmp_edge[0].equals(gene)) {
										from = m2;
									}
									if (tmp_edge[1].equals(gene)) {
										to = m2;
									}
								}
								if (from != -1 && to != -1) {
									remove_node_Edges[kk][0] = from;
									remove_node_Edges[kk][1] = to;
									kk = kk + 1;
								}
							}
							int N_com = countComponents(new_genes_in_vec.size(), remove_node_Edges);
							if (N_com != 1)
								continue;
							/////////////////////////////////////////////////// end of check connection

							M.remove(to_remove_gene);
							//System.out.println("exclude " + to_remove_gene);
							boo = true;
							break;
						}
						if (!boo)
							break; // If no gene could be removed, then break the node check loop.
					} else {
						break; // no gene qualify r_exclude
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////
				
		
			}
			
			
			// finalize the module:
			Object[] modeule_genes_in = M.keySet().toArray();
			Arrays.sort(modeule_genes_in);
			String gg = "";
			for(int k=0; k<modeule_genes_in.length; k++){
				gg = gg + modeule_genes_in[k] + " ";
			}
			
			Vector<String> module_genes = get_vector(modeule_genes_in);
			Set<Integer> module_edge_set = get_edge_set(gene2PPI, PPI, module_genes);
			String edge_str = "";
			for(Integer k:module_edge_set){
				edge_str = edge_str+":"+k;
			}
			//System.out.println("Final module: "+gg);
			double[] tmp_score_1 = get_module_score(nodeHash_gwas, nodeHash_sc, module_genes, module_score);
			
			if(M.size()>=2)modules.put(gg.trim().replaceAll(" ", ":"), seed_gene+"\t"+edge_str.substring(1)+"\t"+tmp_score_1[0]+"\t"+tmp_score_1[1]+"\t"+tmp_score_1[2]);
			
		}
		
		return(modules);
	}
	
	
	public static double[] get_size_statistics (Hashtable<String, Double> nodeHash_gwas, Hashtable<String, Double> nodeHash_sc, Vector<String[]> PPI, Hashtable<String, Vector<Integer>> gene2PPI, int size, String module_score){
		double[] result = new double[2];
		int N = 0;
		
		Object[] genes = nodeHash_gwas.keySet().toArray();
		
		Random r = new Random();
		double[] random_modules_score = new double[1000];
		while(N < 1000) {
			int j = r.nextInt(genes.length);
			String seed_gene = (String)genes[j];
			
			if(!nodeHash_sc.containsKey(seed_gene))continue;
			if(!gene2PPI.containsKey(seed_gene))continue;
			
			Hashtable<String, String> M = new Hashtable<String, String>();
			
			M.put(seed_gene, "");
			while(true){
				Object[] genes_in = M.keySet().toArray();
				if(genes_in.length >= size)break;
				
				Hashtable<String, String> can_nodes = new Hashtable<String, String>();
				for(int n=0; n<genes_in.length; n++){
					String gene = (String)genes_in[n];
					Vector<Integer> cur_idxVec = (Vector<Integer>)gene2PPI.get(gene);
					
					for(Integer k:cur_idxVec) {
						String[] tmp_edge = (String[]) PPI.get(k);
						if(nodeHash_sc.containsKey(tmp_edge[0]) && nodeHash_sc.containsKey(tmp_edge[1])) {
							can_nodes.put(tmp_edge[0], "");
							can_nodes.put(tmp_edge[1], "");
						}
					}
				}
				
				Object[] can_nodes_cur = can_nodes.keySet().toArray();
				int random_select_j = r.nextInt(can_nodes.size());
				M.put( (String)can_nodes_cur[random_select_j], "");
				
			}
			
			
			// finalize the module:
			Object[] modeule_genes_in = M.keySet().toArray();
			double[] tmp_score_1 = get_module_score(nodeHash_gwas, nodeHash_sc, modeule_genes_in, module_score);
			random_modules_score[N] = tmp_score_1[0];
			N++;
		}
		
		result[0] = get_mean(random_modules_score);
		result[1] = get_sd(random_modules_score);
		return(result);
	}
	
	public static double get_mean(double[] arrayD) {
		double mean = 0;
		for(int i=0; i< arrayD.length; i++) {
			mean = mean + arrayD[i];
		}
		mean = mean/arrayD.length;
		return(mean);
	}
	
	public static double get_sd(double[] arrayD) {
		double mean = 0, sd = 0;
		for(int i=0; i< arrayD.length; i++) {
			mean = mean + arrayD[i];
		}
		mean = mean/arrayD.length;
		for(int i=0; i< arrayD.length; i++) {
			sd = sd + Math.pow(arrayD[i]-mean, 2);
		}
		sd = sd/arrayD.length;
		
		return(sd);
	}
	
	public static double[] get_module_score (Hashtable<String, Double> nodeHash_gwas, Hashtable<String, Double> nodeHash_scRNAseq, Vector<String> module_genes, String ms){
		double module_score = 0;
		double denominator = Math.sqrt(module_genes.size());
		
		double z_node_gwas = 0, z_node_scRNAseq = 0;
		for(int i=0; i<module_genes.size(); i++){
			String gene = (String)module_genes.get(i);
			double node_score = 0;
			if(nodeHash_gwas.containsKey(gene))node_score = (double)nodeHash_gwas.get(gene);
			z_node_gwas = z_node_gwas + node_score;
		}
		z_node_gwas = z_node_gwas/denominator;
		
		
		for(int i=0; i<module_genes.size(); i++){
			String gene = (String)module_genes.get(i);
			double node_score = 0;
			if(nodeHash_scRNAseq.containsKey(gene))node_score = (double)nodeHash_scRNAseq.get(gene);
			z_node_scRNAseq = z_node_scRNAseq + node_score;
		}
		z_node_scRNAseq = z_node_scRNAseq/denominator;
		
		if(ms.equals("tw")) {
			module_score = z_node_gwas + z_node_scRNAseq;
		} else if(ms.equals("penalty")) {
			double m = (z_node_gwas + z_node_scRNAseq)/2;
			module_score = z_node_gwas + z_node_scRNAseq - Math.sqrt( (z_node_gwas - m)*(z_node_gwas - m) + (z_node_scRNAseq - m)*(z_node_scRNAseq - m) );
		}
		
		double[] res = new double[3];
		res[0] = module_score;
		res[1] = z_node_gwas;
		res[2] = z_node_scRNAseq;
		
		return(res);
	}
	
	public static double[] get_module_score (Hashtable<String, Double> nodeHash_gwas, Hashtable<String, Double> nodeHash_scRNAseq, Object[] module_genes, String ms){
		double module_score = 0;
		double denominator = Math.sqrt(module_genes.length);
		
		double z_node_gwas = 0, z_node_scRNAseq = 0;
		for(int i=0; i<module_genes.length; i++){
			String gene = (String)module_genes[i];
			double node_score = 0;
			if(nodeHash_gwas.containsKey(gene))node_score = (double)nodeHash_gwas.get(gene);
			z_node_gwas = z_node_gwas + node_score;
		}
		z_node_gwas = z_node_gwas/denominator;
		
		
		for(int i=0; i<module_genes.length; i++){
			String gene = (String)module_genes[i];
			double node_score = 0;
			if(nodeHash_scRNAseq.containsKey(gene))node_score = (double)nodeHash_scRNAseq.get(gene);
			z_node_scRNAseq = z_node_scRNAseq + node_score;
		}
		z_node_scRNAseq = z_node_scRNAseq/denominator;
		
		if(ms.equals("tw")) {
			module_score = z_node_gwas + z_node_scRNAseq;
		} else if(ms.equals("penalty")) {
			double m = (z_node_gwas + z_node_scRNAseq)/2;
			module_score = z_node_gwas + z_node_scRNAseq - Math.sqrt( (z_node_gwas - m)*(z_node_gwas - m) + (z_node_scRNAseq - m)*(z_node_scRNAseq - m) );
		}
		
		double[] res = new double[3];
		res[0] = module_score;
		res[1] = z_node_gwas;
		res[2] = z_node_scRNAseq;
		
		return(res);
	}
	
	public static Vector<Hashtable<String, Double>> get_normalization_calibration (
			Hashtable<String, Double> nodeHash_gwas, 
			Hashtable<String, Double> nodeHash_sc, 
			boolean verbose, 
			String outFile, 
			String cell_type, int cell_type_index) throws IOException{
		
		Vector<String> shared_genes   = new Vector<String>();
		Hashtable<String, Double> share_nodeHash_gwas = new Hashtable<String, Double>();
		Hashtable<String, Double> share_nodeHash_scrn = new Hashtable<String, Double>();
		
		double gwas_sd = 0, gwas_mean = 0;
		Object[] genes = nodeHash_gwas.keySet().toArray();
		for(int j=0; j<genes.length; j++){
			double node_score = (double)nodeHash_gwas.get(genes[j]);
			gwas_mean = gwas_mean + node_score;
			
			if(nodeHash_sc.containsKey(genes[j]))shared_genes.add((String)genes[j]);
		}
		gwas_mean = gwas_mean/genes.length;
		for(int j=0; j<genes.length; j++){
			double node_score = (double)nodeHash_gwas.get(genes[j]);
			gwas_sd = gwas_sd + Math.pow(node_score - gwas_mean, 2);
		}
		gwas_sd = Math.sqrt(gwas_sd/(genes.length-1) );
		System.out.println("Raw GWAS genes: n = "+nodeHash_gwas.size()+", mean: "+String.format("%6.4e",gwas_mean)+", sd = "+String.format("%6.4e",gwas_sd));
		
		
		double scrn_sd = 0, scrn_mean = 0;
		genes = nodeHash_sc.keySet().toArray();
		for(int j=0; j<genes.length; j++){
			double node_score = (double)nodeHash_sc.get(genes[j]);
			scrn_mean = scrn_mean + node_score;
		}
		scrn_mean = scrn_mean/genes.length;
		for(int j=0; j<genes.length; j++){
			double node_score = (double)nodeHash_sc.get(genes[j]);
			scrn_sd = scrn_sd + Math.pow(node_score - scrn_mean, 2);
		}
		scrn_sd = Math.sqrt(scrn_sd/(genes.length-1) );
		System.out.println("Raw scRNAseq genes: n = "+nodeHash_sc.size()+", mean: "+String.format("%6.4e",scrn_mean)+", sd = "+String.format("%6.4e",scrn_sd));
		
		// get mean and sd for shared genes
		gwas_sd = 0; gwas_mean = 0; scrn_sd = 0; scrn_mean = 0;
		for(int j=0; j<shared_genes.size(); j++){
			String gene = shared_genes.get(j);
			gwas_mean = gwas_mean + (double)nodeHash_gwas.get(gene);
			scrn_mean = scrn_mean + (double)nodeHash_sc.get(gene);
		}
		scrn_mean = scrn_mean/shared_genes.size();
		gwas_mean = gwas_mean/shared_genes.size();
		
		for(int j=0; j<shared_genes.size(); j++){
			String gene = shared_genes.get(j);
			gwas_sd = gwas_sd + Math.pow((double)nodeHash_gwas.get(gene) - gwas_mean, 2);
			scrn_sd = scrn_sd + Math.pow((double)nodeHash_sc.get(gene) - scrn_mean, 2);
		}
		scrn_sd = Math.sqrt(scrn_sd/(shared_genes.size()-1) );
		gwas_sd = Math.sqrt(gwas_sd/(shared_genes.size()-1) );
		System.out.println("Shared genes: n = "+shared_genes.size());
		System.out.println("Shared genes in GWAS mean: "+String.format("%6.4e",gwas_mean)+", sd = "+String.format("%6.4e",gwas_sd));
		System.out.println("Shared genes in scRNAseq mean: "+String.format("%6.4e",scrn_mean)+", sd = "+String.format("%6.4e",scrn_sd));
		System.out.println("Generate standard normal distribution for "+shared_genes.size()+" genes");
		
		double[] ref_gwas = QuantileNormalization.generate_standard_normal_distribution(shared_genes.size());
		double[][] raw_data = new double[shared_genes.size()][3];
		for(int j=0; j<shared_genes.size(); j++){
			String gene = shared_genes.get(j);
			double gwas_score = ((double)nodeHash_gwas.get(gene) - gwas_mean)/gwas_sd;
			double scrn_score = ((double)nodeHash_sc.get(gene) - scrn_mean)/scrn_sd;
			raw_data[j][0] = gwas_score;
			raw_data[j][1] = scrn_score;
			
			raw_data[j][2] = ref_gwas[j];
		}
		double[][] raw_data_clone = raw_data.clone();
		
		// quantile normalization
		double[][] quantile_normalized_data = QuantileNormalization.quantilenormalize(raw_data);
		
		
		// shift to zero-centered
		scrn_mean = 0; gwas_mean = 0;
		for(int i=0; i<shared_genes.size(); i++) {
			gwas_mean = gwas_mean + quantile_normalized_data[i][0];
			scrn_mean = scrn_mean + quantile_normalized_data[i][1];
		}
		scrn_mean = scrn_mean/shared_genes.size();
		gwas_mean = gwas_mean/shared_genes.size();
				
				
		// assign the quantile-normalized data into hashtables
		for(int j=0; j<shared_genes.size(); j++){
			String gene = shared_genes.get(j);
			share_nodeHash_gwas.put(gene, quantile_normalized_data[j][0]-gwas_mean);
			share_nodeHash_scrn.put(gene, quantile_normalized_data[j][1]-scrn_mean);
		}
		
		
		if(verbose) {
			File f = new File(outFile+"..shared_genes.scores.calibration.txt");
			if(cell_type_index==0 & f.exists()) {
				f.delete();
			}
			
			if(cell_type_index==0) {
				BufferedWriter bout = new BufferedWriter(new FileWriter(f));
				bout.write("gene\tgwas_score_raw\tscrn_score_raw\tref_random_score\tgwas_score_raw_z\tscrn_score_raw_z\tref_calibrated\tgwas_score_raw_z_calibrated\tscrn_score_raw_z_calibrated\tcell_type\n");
				for(int j=0; j<shared_genes.size(); j++){
					String gene = shared_genes.get(j);
					bout.write(gene+"\t"+nodeHash_gwas.get(gene)+"\t"+nodeHash_sc.get(gene)+"\t"+raw_data[j][2]+"\t"+raw_data_clone[j][0]+"\t"+raw_data_clone[j][1]+"\t"+quantile_normalized_data[j][2]+"\t"+
					          share_nodeHash_gwas.get(gene)+"\t"+share_nodeHash_scrn.get(gene)+"\t"+cell_type+"\n");
					
				} bout.flush(); bout.close();
			} else {
				BufferedWriter bout = new BufferedWriter(new FileWriter(f, true));
				for(int j=0; j<shared_genes.size(); j++){
					String gene = shared_genes.get(j);
					bout.write(gene+"\t"+nodeHash_gwas.get(gene)+"\t"+nodeHash_sc.get(gene)+"\t"+raw_data[j][2]+"\t"+raw_data_clone[j][0]+"\t"+raw_data_clone[j][1]+"\t"+quantile_normalized_data[j][2]+"\t"+
					          share_nodeHash_gwas.get(gene)+"\t"+share_nodeHash_scrn.get(gene)+"\t"+cell_type+"\n");
					
				} bout.flush(); bout.close();
			}
			
		}
		
		Vector<Hashtable<String, Double>> vec = new Vector<Hashtable<String, Double>>();
		vec.add(share_nodeHash_gwas);
		vec.add(share_nodeHash_scrn);
		return(vec);
	}

	
	public static Hashtable<String, Double> get_rank_score_single (Hashtable<String, Double> nodeHash_gwas){
		Hashtable<String, Double> rank_nodeHash_gwas = new Hashtable<String, Double>();
		Object[] genes = nodeHash_gwas.keySet().toArray();
		int N = nodeHash_gwas.size();
		double[] gwas_score = new double[N];
		
		for(int j=0; j<genes.length; j++){
			gwas_score[j] = (double)nodeHash_gwas.get(genes[j]);
		}
		
		double[] gwas_rank = rankify.rankify(gwas_score);
		
		for(int j=0; j<genes.length; j++){
			String gene = (String)genes[j];
			rank_nodeHash_gwas.put(gene, gwas_rank[j]/N);
		}
		
		return(rank_nodeHash_gwas);
	}

	public static Vector<String> get_vector(Object[] input){
		Vector<String> output = new Vector<String>();
		for(int i=0; i<input.length; i++) {
			output.add((String)input[i]);
		}
		return(output);
	}
	
	
	public static Set<Integer> get_edge_set(Hashtable<String, Vector<Integer>> gene2PPI, Vector<String[]> PPI, Vector<String> nodes) {
		Set<Integer> module_edge_set = new HashSet<Integer>(new Vector<Integer>());
		
		for(int n=0; n<nodes.size(); n++){
			String tmp_gene = (String)nodes.get(n);
			if(gene2PPI.containsKey(tmp_gene)){
				Vector<Integer> cur_idxVec = (Vector<Integer>)gene2PPI.get(tmp_gene);
				
				for(int n2=0; n2 < cur_idxVec.size(); n2++){
					int N2 = (Integer)cur_idxVec.get(n2);
					String[] tmp_edge_tmp = (String[]) PPI.get(N2);
					
					if(is_element(tmp_edge_tmp[0], nodes) && is_element(tmp_edge_tmp[1], nodes) ){
						module_edge_set.add(N2);
					}
				}
			}
		}
		
		return(module_edge_set);
	}
	
	
	public static boolean is_element(String a, Vector<String> b) {
		boolean check = false;
		for(int i=0; i<b.size(); i++){
			String str = (String)b.get(i);
			if(a.equals(str)){
				check = true;
				break;
			}
		}
		return(check);
	}
	
	public static int countComponents(int n, int[][] edges) {
		int count = n;

		int[] root = new int[n];
		
		for (int i = 0; i < n; i++) {
			root[i] = i;
		}

		for (int i = 0; i < edges.length; i++) {
			int x = edges[i][0];
			int y = edges[i][1];

			int xRoot = getRoot(root, x);
			int yRoot = getRoot(root, y);

			if (xRoot != yRoot) {
				count--;
				root[xRoot] = yRoot;
			}
		}
		return count;
	}
	
	public static int getRoot(int[] arr, int i) {
		while (arr[i] != i) {
			arr[i] = arr[arr[i]];
			i = arr[i];
		}
		return i;
	}

}
