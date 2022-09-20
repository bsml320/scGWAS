package util;

import java.io.IOException;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;

public class r_to_z {
	
	/*
	public static void main() throws IOException {
		//chisq2z(new double[] {1.64});
	}
	*/
	public static double[] transform (double[] pcc, double NF) throws IOException{
		double pcc_to_t[] = new double[pcc.length];
		for(int i=0; i<pcc.length; i++) {
			pcc_to_t[i] = pcc[i] * Math.sqrt(NF-2)/Math.sqrt(1-pcc[i] * pcc[i]);
		}
		
		double mean_pcc_to_t = 0, sd_pcc_to_t = 0;
		for(int i=0; i<pcc.length; i++) {
			mean_pcc_to_t = mean_pcc_to_t + pcc_to_t[i];
		}
		mean_pcc_to_t = mean_pcc_to_t/pcc_to_t.length;
		System.out.println("mean_pcc_to_t = "+mean_pcc_to_t);
		
		for(int i=0; i<pcc.length; i++) {
			sd_pcc_to_t = sd_pcc_to_t + Math.pow(pcc_to_t[i] - mean_pcc_to_t, 2);
		}
		sd_pcc_to_t = Math.sqrt(sd_pcc_to_t/(pcc_to_t.length-1) );
		System.out.println("sd_pcc_to_t = "+sd_pcc_to_t);
		
		double degf = NF-2;
		double adjust_sigma =  sd_pcc_to_t * Math.sqrt(degf-2)/Math.sqrt(degf); 
		double normalized_PCC_to_t[] = new double[pcc.length];
		for(int i=0; i<pcc.length; i++) {
			normalized_PCC_to_t[i] = Math.abs(pcc_to_t[i]-mean_pcc_to_t)/adjust_sigma;
		}
		System.out.println("adjust_sigma = "+adjust_sigma);
		
		TDistribution tDistribution = new TDistribution(degf);
		double p_from_t[] =  new double[pcc.length];
		for(int i=0; i<pcc.length; i++) {
			p_from_t[i] = 2 * (1 - tDistribution.cumulativeProbability( normalized_PCC_to_t[i] ));
		}
		
		double z_[] =  new double[pcc.length];
		double max = 0;
		NormalDistribution d = new NormalDistribution (); 
		for(int i=0; i<pcc.length; i++) {
			z_[i] = d.inverseCumulativeProbability(1-p_from_t[i]);
			if(!Double.isInfinite(z_[i])) {
				if(z_[i] > max)max = z_[i];
			}
		}
		
		for(int i=0; i<pcc.length; i++) {
			if(Double.isInfinite(z_[i])) {
				z_[i] = max ;
			}
		}
		
		return(z_);
	}
	
	
	public static double[] Fisher (double[] pcc) throws IOException{
		double pcc_to_z[] = new double[pcc.length];
		for(int i=0; i<pcc.length; i++) {
			pcc_to_z[i] = 0.5 * ( Math.log(1+pcc[i]) - Math.log(1-pcc[i]) );
		}
		
		return(pcc_to_z);
	}
	
	
	public static double[] chisq2z(double[] chisq) throws IOException{
		ChiSquaredDistribution cd = new ChiSquaredDistribution(1);
		NormalDistribution d = new NormalDistribution (); 
		double[] z_ = new double[chisq.length];
		double max = 0;
		for(int i=0; i<chisq.length; i++) {
			double p = cd.cumulativeProbability(chisq[i]);
			double z = d.inverseCumulativeProbability(p);
			z_[i] = z;
			if(!Double.isInfinite(z_[i])) {
				if(z_[i] > max)max = z_[i];
			}
			System.out.println(chisq[i]+"\t"+p+"\t"+z);
		}
		
		for(int i=0; i<chisq.length; i++) {
			if(Double.isInfinite(z_[i])) {
				z_[i] = max ;
			}
		}
		return(z_);
	}
	
}
