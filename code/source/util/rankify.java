package util;

public class rankify {
	public static double[] rankify(double A[]) {
		int N = A.length;
		double[] Rank_X = new double[N];
		
		for(int i=0; i<N; i++) {
			int r = 1, s = 1;
			for(int j=0; j<i; j++) {
				if(A[j] < A[i]) {
					r = r + 1;
				}
				if(A[j] == A[i]) {
					s = s + 1;
				}
			}
			
			for(int j=i+1; j<N; j++) {
				if(A[j] < A[i]) {
					r = r + 1;
				}
				if(A[j] == A[i]) {
					s = s + 1;
				}
			}
			Rank_X[i] = r + (s-1) * 0.5;
		}
		return Rank_X;
	}
	/*
	public static double[] rankify(double A[], int n) {
		// Rank Vector
		double R[] = new double[n];
		
		for (int i = 0; i < n; i++) {
			int r = 1, s = 1;

			for (int j = 0; j < n; j++) {
				if (j != i && A[j] < A[i])
					r += 1;

				if (j != i && A[j] == A[i])
					s += 1;
			}

			// obtain rank
			R[i] = r + (float) (s - 1) / (float) 2;
		}
		return(R);
	}
	*/
}
