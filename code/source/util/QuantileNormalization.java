package util;

import java.util.*;

import org.apache.commons.math3.distribution.NormalDistribution;

public class QuantileNormalization {
	public static double[][] quantilenormalize(double[][] rawData) {
       
       
        //Calculate the average expression, when per sample all raw expression levels have been ordered:
        int probeCount = rawData.length;
        int sampleCount = rawData[probeCount-1].length;
        System.out.println("Performing quantile normalization for "+probeCount + " rows and " +sampleCount+" columns");
        double[][] qn_data = new double[probeCount][sampleCount];
        
        double[] rankedMean = new double[probeCount];
        for (int sampleID=0; sampleID<sampleCount; sampleID++) {
            double[] x = new double[probeCount]; 
            
            for (int probeID=0; probeID<probeCount; probeID++) {
                x[probeID] = rawData[probeID][sampleID];
            }
            java.util.Arrays.sort(x);
            
            for (int probeID=0; probeID<probeCount; probeID++) {
                rankedMean[probeID] += x[probeID];
            }
        }
        
        for (int probeID=0; probeID<probeCount; probeID++) {
            rankedMean[probeID]/=(double) sampleCount;
        }
        
        //Iterate through each sample: skip s=0 because this is the column for reference. Won't be used any way.
        for (int s=0; s<sampleCount; s++) {
            double[] probes = new double[probeCount]; 
            for (int p=0; p<probeCount; p++) {
                probes[p]=rawData[p][s];
            }
            
            double[] probesRanked = rankify.rankify(probes);
            
            double[] probesQuantileNormalized = new double[probeCount];
            for (int p=0; p<probeCount; p++) {
                probesQuantileNormalized[p] = rankedMean[ (int)probesRanked[p] - 1 ];
            }
            for (int p=0; p<probeCount; p++) {
            	qn_data[p][s] = (float) probesQuantileNormalized[p];
            }
        }
        
        return(qn_data);
    }
	
	public static double[][] quantilenormalize_2cn(double[] refData, double[] rawData) {

		// Calculate the average expression, when per sample all raw expression levels
		// have been ordered:
		int probeCount = rawData.length;
		int sampleCount = 2;
		System.out.println(
				"Performing quantile normalization for " + probeCount + " rows and " + sampleCount + " columns");
		double[][] qn_data = new double[probeCount][sampleCount];

		double[] rankedMean = new double[probeCount];

		java.util.Arrays.sort(refData);
		java.util.Arrays.sort(rawData);
		for (int probeID = 0; probeID < probeCount; probeID++) {
			rankedMean[probeID] = rankedMean[probeID] + refData[probeID] + rawData[probeID];
		}

		for (int probeID = 0; probeID < probeCount; probeID++) {
			rankedMean[probeID] /= (double) sampleCount;
		}

		// Iterate through each sample: skip s=0 because this is the column for
		// reference. Won't be used any way.
		double[] probes = new double[probeCount];
		for (int p = 0; p < probeCount; p++) {
			probes[p] = rawData[p];
		}

		double[] probesRanked = rankify.rankify(probes);

		double[] probesQuantileNormalized = new double[probeCount];
		for (int p = 0; p < probeCount; p++) {
			probesQuantileNormalized[p] = rankedMean[(int) probesRanked[p] - 1];
		}
		for (int p = 0; p < probeCount; p++) {
			qn_data[p][1] = (float) probesQuantileNormalized[p];
		}

        return(qn_data);
    }
	
	public static int[] rank_double_array(double[] original) {
		int N = original.length;
		int[] ranked = new int[N];

		// create an empty TreeMap
		Map<Double, Integer> map = new TreeMap<>();

		// store (element, index) pair in TreeMap
		for (int i = 0; i < N; i++) {
			map.put(original[i], i);
		}
		// keys are stored in sorted order in TreeMap

		// rank starts from 1
		int rank = 0;

		// iterate through the map and replace each element by its rank
		
		for (Map.Entry<Double, Integer> entry : map.entrySet()) {
			ranked[entry.getValue()] = rank++;
		}
		System.out.println("total rank: "+rank);
		return(ranked);
	}
	
	public static int[] rank_double_array_2(double[] original) {
		int N = original.length;
		int[] ranked = new int[N];

		for(int i=0; i<N; i++) {
			int count = 0;
			for(int j=0; j<N; j++) {
				if(j==i)continue;
				if(original[j] < original[i])count++;
			}
			ranked[i] = count;
		}
		
		return(ranked);
	}
	
	
	public static double[] generate_normal_distribution(int N) {
		double[] norm = new double[N];
		Random r = new Random();
		
		for(int i=0; i<N; i++) {
			norm[i] = r.nextGaussian();
		}
		
		double sd = 0, mean = 0;
		for(int i=0; i<N; i++) {
			mean = mean + norm[i];
		}
		mean = mean/N;
		for(int i=0; i<N; i++) {
			sd = sd + Math.pow(norm[i] - mean, 2);
		}
		sd = Math.sqrt(sd/(N-1) );
		System.out.println("Random scores: n = "+N+", mean: "+String.format("%6.4e",mean)+", sd = "+String.format("%6.4e",sd));
		
		return(norm);
	}
    
	public static double[] generate_standard_normal_distribution(int N) {
		double[] norm = new double[N];
		
		NormalDistribution d = new NormalDistribution (); 
		for(int i=0; i<N; i++) {
			double di = i;
			double p_rank = di/N;
			if(i == 0)p_rank = 0.5/N;
			norm[i] = d.inverseCumulativeProbability(p_rank);
			if(i==0)System.out.println(N+"\t"+p_rank+"\t"+norm[i]);
		}
		
		//System.out.println(d.inverseCumulativeProbability(0.05));
		
		double sd = 0, mean = 0;
		for(int i=0; i<N; i++) {
			mean = mean + norm[i];
		}
		mean = mean/N;
		for(int i=0; i<N; i++) {
			sd = sd + Math.pow(norm[i] - mean, 2);
		}
		sd = Math.sqrt(sd/(N-1) );
		System.out.println("Random scores: n = "+N+", mean: "+String.format("%6.4e",mean)+", sd = "+String.format("%6.4e",sd));
		
		return(norm);
	}
}
