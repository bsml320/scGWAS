package main;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public class scGWAS {
	public static void main(String[] args) {
		//System.out.println("Hello World");
		try {
			String  gwas_node_file      = null;
			String  scrn_expr_file      = null;
			String  network_file        = null;
			String  network_weight_file = null;
			String  outfile             = null;
			String  run_model           = "node";
			String  module_score        = "penalty";
			String  normalization_model = "calibration";
			double  r_include           = 0.2;
			double  r_exclude           = 0.1;
			boolean permutation         = true;
			boolean remove_zero         = false;
			boolean verbose             = false;
			String  cell_type_str       = null;
			Hashtable<String, String> cell_type_ht = new Hashtable<String, String>();
			String int_genes_file       = null;
			String exclude_genes_file   = null;
			int module_num              = 1000;
			
			Hashtable<String, String> inputHT = null;
			
			if(args.length==1) {
				File f = new File(args[0]);
				if(f.exists()) {
					System.out.println("Read parameters from the configure file: "+args[0]);
					inputHT = new Hashtable<String, String>();
					BufferedReader bin = new BufferedReader (new FileReader (f));
					String str;
					while((str=bin.readLine())!=null) {
						if(str.startsWith("#"))continue;
						String[] tmp = str.split("=");
						if(tmp.length>=2) {
							inputHT.put(tmp[0], tmp[1]);
						}
					}bin.close();
				} else if(args[0].equals("--help")){
					System.out.println("usage: java -jar scGWAS_v1.jar configure.txt");
					System.out.println("run_model=[node|one_node|reference]\tRequired. This is the task you want to perform.");
					System.out.println("gwas_node_file=file_location\tRequired if run_model=node or node_edge. This is the a 2-column, tab-seperated file with gene-based z-scores from GWAS. The first column is gene name. The second column is z-score.");
					System.out.println("scrn_expr_file=file_location\tRequired. This is the scRNAseq expression file. Genes on rows and cells (or cell types) on columns. Header required.");
					System.out.println("network_file=file_location\tRequired if run_model=node|one_node|reference. PPI file. Each row includes one pair of genes. Tab-seperated. Gene name required (not gene ID).");
					System.out.println("outfile=file_location.\tOptional but highly suggested to provide. If not specified, the output parameter will be automatically assigned using the current time.");
					System.out.println("normalization_model=[scale|rank|calibration]\tRequired if run_model=node|one_node|reference. The normalization method.");
					System.out.println("r_include=0.1\tOptional. The threshold to include a neighbor node. Defalut 0.1");
					System.out.println("r_exclude=0.05\tOptional. The threshold to exclude a component gene. Defalut 0.05");
					System.out.println("permutation=[false|true]\tOptional. Case sensitive.");
					System.out.println("verbose=[true|false]\tOptional.");
					System.out.println("cell_type_str=[ct1,ct2]\tOptional. If assigned, scGWAS will only consider the cell types listed in this string. Use comma to seperate cell types with no space.");
					
					System.exit(0);
				} else {
					System.out.println("File "+args[0]+" does not exist.");
					System.exit(0);
				}
			} else {
				// 2. input strings
				inputHT = getinput(args);
			}
			
			
			if(!inputHT.containsKey("run_model")) {
				System.out.println("Must provide run_model: node, one_node, reference.");
				System.exit(0);
			}
			
			Hashtable<String, String> effective_model = new Hashtable<String, String>();
			effective_model.put("node", "");
			effective_model.put("one_node", "");
			effective_model.put("reference", "");
			
			run_model = (String)inputHT.get("run_model");
			System.out.println("run_model="+run_model);
			if(!effective_model.containsKey(run_model)) {
				System.out.println("run_model must be one of <node, one_node, reference>");
				System.exit(0);
			}
			
			System.out.println("The following parameters will be used: ");
			if(inputHT.containsKey("gwas_node_file")){
				gwas_node_file = (String)inputHT.get("gwas_node_file");
				System.out.println("gwas_node_file="+gwas_node_file);
			}
			
			if(inputHT.containsKey("scrn_expr_file")){
				scrn_expr_file = (String)inputHT.get("scrn_expr_file");
				System.out.println("scrn_expr_file="+scrn_expr_file);
			}
			
			if(inputHT.containsKey("module_score")){
				module_score = (String)inputHT.get("module_score");
			} else {
				System.out.println("module_score set as tw for <two node weights>.");
			}
			
			if(inputHT.containsKey("network_file")){
				network_file = (String)inputHT.get("network_file");
				System.out.println("network_file="+network_file);
			}
			
			if(inputHT.containsKey("network_weight_file")){
				network_weight_file = (String)inputHT.get("network_weight_file");
			}
			
			if(inputHT.containsKey("outfile")){
				outfile =(String)inputHT.get("outfile");
				System.out.println("outfile="+outfile);
			} else {
				Date dNow = new Date();
				DateFormat dateFormat = new SimpleDateFormat("yyyy-mm-dd-hh-mm-ss");
				outfile = dateFormat.format(dNow);
				System.out.println("Set outfile="+outfile);
			}
			
			if(inputHT.containsKey("normalization_model")){
				normalization_model = (String)inputHT.get("normalization_model");
				System.out.println("normalization_model="+normalization_model);
			}
			
			if(inputHT.containsKey("r_include")){
				r_include = Double.parseDouble((String)inputHT.get("r_include"));
				System.out.println("r_include="+r_include);
			}
			
			if(inputHT.containsKey("r_exclude")){
				r_exclude = Double.parseDouble((String)inputHT.get("r_exclude"));
				System.out.println("r_exclude="+r_exclude);
			}
			
			if(inputHT.containsKey("verbose")){
				verbose = java.lang.Boolean.parseBoolean((String)inputHT.get("verbose"));
			}
			
			if(inputHT.containsKey("cell_type_str")){
				cell_type_str = (String)inputHT.get("cell_type_str");
				String[] line = cell_type_str.split(",");
				for(int i=0; i<line.length; i++)cell_type_ht.put(line[i], "");
				System.out.println("Read in "+line.length+" cell types.");
			} else {
				cell_type_ht = null;
				if(run_model.equalsIgnoreCase("node") | run_model.equalsIgnoreCase("node_edge")) System.out.println("No cell type specified. Using all cell types in the input file.");
			}
			
			if(inputHT.containsKey("permutation")){
				permutation = java.lang.Boolean.parseBoolean((String)inputHT.get("permutation"));
				System.out.println("permutation="+permutation);
				if(!permutation) {
					System.out.println("No permutation will be performed");
				}
				
			} else {
				permutation = false;
				System.out.println("No permutation will be performed");
			}
			
			if(inputHT.containsKey("remove_zero")){
				remove_zero = java.lang.Boolean.parseBoolean((String)inputHT.get("remove_zero"));
			}
			
			if(inputHT.containsKey("int_genes_file")){
				int_genes_file = inputHT.get("int_genes_file");
				System.out.println("int_genes_file="+int_genes_file);
			}
			
			if(inputHT.containsKey("exclude_genes_file")){
				exclude_genes_file = inputHT.get("exclude_genes_file");
				System.out.println("exclude_genes_file="+exclude_genes_file);
			}
			
			if(inputHT.containsKey("module_num")){
				module_num = Integer.parseInt(inputHT.get("module_num"));
				System.out.println("module_num="+module_num);
			}
			/////////////////
			
			
			
			if(run_model.equalsIgnoreCase("node")) {
				System.out.println("Setting run_model="+run_model+". Two-way node-weighted dense module construction.");
				
				Date dNow = new Date();
				SimpleDateFormat ft = new SimpleDateFormat("E yyyy.MM.dd 'at' hh:mm:ss a zzz");
				System.out.println("Start at: " + ft.format(dNow));
				
				if(scrn_expr_file==null | gwas_node_file==null | network_file==null) {
					System.out.println("When run_model=node, the following parameters are required: gwas_node_file, scrn_expr_file, and network_file.");
					System.exit(0);
				}
				
				if(network_weight_file!=null) {
					System.out.println("network_weight_file="+network_weight_file);
					two_nodes_scGWAS_weight.scGWAS(gwas_node_file, scrn_expr_file, module_score, network_file, network_weight_file, normalization_model, cell_type_ht, r_include, r_exclude, outfile, int_genes_file, exclude_genes_file, permutation, remove_zero, module_num, verbose);
				} else {
					two_nodes_scGWAS.scGWAS_wrapper(gwas_node_file, scrn_expr_file, module_score, network_file, normalization_model, cell_type_ht, r_include, r_exclude, outfile, int_genes_file, exclude_genes_file, permutation, remove_zero, module_num, verbose);
				}
				
				System.out.println("Finished: " + ft.format(new Date()));
			}
			
			// Tag #2 one node
			if(run_model.equalsIgnoreCase("one_node")) {
				System.out.println("Setting run_model="+run_model+". One-way node-weighted dense module construction.");
				
				Date dNow = new Date();
				SimpleDateFormat ft = new SimpleDateFormat("E yyyy.MM.dd 'at' hh:mm:ss a zzz");
				System.out.println("Start at: " + ft.format(dNow));
				
				one_node_MEBE.one_node_(scrn_expr_file, network_file, cell_type_ht, r_include, r_exclude, outfile, int_genes_file, exclude_genes_file, permutation, module_num, verbose);
				
				System.out.println("Finished: " + ft.format(new Date()));
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static Hashtable<String, String> getinput (String[] args) throws Exception{
		Hashtable<String, String> inputHt = new Hashtable<String, String>();
		for(int i=0; i<args.length; i++){
			if(args[i].startsWith("#"))continue;
			if(args[i].indexOf("=")==-1){
				System.out.println("The following input parameters have problems: "+args[i]);
				System.exit(0);
			}
			String[] tmp = args[i].split("=");
			if(tmp.length!=2){
				System.out.println("The following input parameters have problems: "+args[i]);
				System.exit(0);
			}
			inputHt.put(tmp[0], tmp[1]);
		}
		
		return(inputHt);
	}
	
}
