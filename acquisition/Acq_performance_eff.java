
package gps.acquisition;

import cgra.pe.PETrigonometry;

public class Acquisition {

	// Static input parameters
	private static final int FREQ_STEP_COUNT = 11;
	private static final int FREQ_MOD_MIN = -5000;
	private static final int FREQ_MOD_MAX = 5000;
	private static final int FREQ_MOD_STEP = 1000;
	private static final int FREQ_SAMPLING = 400000;
	private static final float ACQ_THRESHOLD = 0.015f;
	private static final float FPI = (float) 3.141;

	// Results
	private int dopplerShift = 0;
	private int codeShift = 0;

	// Dynamic input parameters
	private final int sCount;

	// Sample storage for unmodified samples
	private final float[] sample_r;
	private final float[] sample_i;
	private int sampleWriteIndex = 0;

	// Storage for intermediate results
	private final float[] res_r_1;
	private final float[] res_i_1;
	private final float[] res_r_2;
	private final float[] res_i_2;

	// Code storage (not transformed)
	private final float[] code_r;
	private final float[] code_i;
	private int codeWriteIndex = 0;

	public Acquisition(int nrOfSamples) {
		this.sCount = nrOfSamples;

		this.sample_r = new float[this.sCount];
		this.sample_i = new float[this.sCount];
		this.res_r_1 = new float[this.sCount];
		this.res_i_1 = new float[this.sCount];
		this.res_r_2 = new float[this.sCount];
		this.res_i_2 = new float[this.sCount];

		this.code_r = new float[this.sCount];
		this.code_i = new float[this.sCount];
	}
	
	public void enterSample(float real, float imag) {
		sample_r[sampleWriteIndex] = real;
		sample_i[sampleWriteIndex] = imag;
		++sampleWriteIndex;
	}
	
	public void enterCode(float real, float imag) {
		code_r[codeWriteIndex] = real;
		code_i[codeWriteIndex] = imag;
		++codeWriteIndex;
	}
	
	public boolean startAcquisition() {
		final float lFPI = FPI;
		final int sCount = this.sCount;
		float sumreal = 0;
		float sumimag = 0;
		
		float pin = 0;
		float smax = 0;
		int fdAtMax = 0;
		int tauAtMax = 0;
		
		// --- Calculate Smax on the fly while performing all the other operations
		for (int j = 0; j < FREQ_STEP_COUNT; j=j+1)
		{
		
			// --- For each fd, compute Xin[n]*e^(-2j*pi*fd*n/fs)
			for (int i = 0; i < sCount; i=i+1)
			{
				res_r_1[i] =  sample_r[i] * PETrigonometry.cos(2 * lFPI/FREQ_SAMPLING * (FREQ_MOD_MIN + FREQ_MOD_STEP * j) * i ) + sample_i[i] * PETrigonometry.sin(2 * lFPI/FREQ_SAMPLING * (FREQ_MOD_MIN + FREQ_MOD_STEP * j) * i );
				res_i_1[i] = -sample_r[i] * PETrigonometry.sin(2 * lFPI/FREQ_SAMPLING * (FREQ_MOD_MIN + FREQ_MOD_STEP * j) * i ) + sample_i[i] * PETrigonometry.cos(2 * lFPI/FREQ_SAMPLING * (FREQ_MOD_MIN + FREQ_MOD_STEP * j) * i );
			}
			
			// circular convolusion. The formula x1 convolution x2 <-> X1*X2 is used
			// the convolution in time domain corresponds to a sum of res_r_1[t]*code_r([k-t]modulo sCount)
			// k-t modulo sCount can also be a negative value. Instead of using conditional clause, the inner loop is divided into 2 innner loops
			// because in frequency domain we have DFT*(code), we need to make an adjustment in time domain (complex conjug. & code[sCount-index] modulo sCount 
			for (int k=0;k<sCount;k=k+1) {
				sumreal = 0;
				sumimag = 0;
				for (int t = 0; t < k; t=t+1)
				{	
					sumreal +=   res_r_1[t]* code_r[(sCount -k+t)] + res_i_1[t]* code_i[(sCount -k+t)]; 
					sumimag +=  -res_r_1[t]* code_i[(sCount -k+t)] + res_i_1[t]* code_r[(sCount -k+t)];
				}
				sumreal +=   res_r_1[k]* code_r[0] + res_i_1[k]* code_i[0]; 
				sumimag +=  -res_r_1[k]* code_i[0] + res_i_1[k]* code_r[0];
				
				for (int t = k+1; t < sCount; t=t+1)
				{	
					sumreal +=   res_r_1[t]* code_r[t-k] + res_i_1[t]* code_i[t-k]; 
					sumimag +=  -res_r_1[t]* code_i[t-k] + res_i_1[t]* code_r[t-k];
				}
				res_r_2[k] 	=  sumreal;
				res_i_2[k] 	=  sumimag;
			}

			// Compute the maximum value
			for (int i = 0; i < sCount; i=i+1) {
				if (smax < res_r_2[i] * res_r_2[i] + res_i_2[i] * res_i_2[i]) {
					smax = (res_r_2[i] * res_r_2[i] + res_i_2[i] * res_i_2[i]);
					fdAtMax = (FREQ_MOD_MIN + FREQ_MOD_STEP * j);
					tauAtMax = i;
				}
			}
		}
		
		// Compute results
		smax = smax/sCount;
		dopplerShift =  fdAtMax;
		codeShift = tauAtMax;
		return (smax*sCount / pin) > ACQ_THRESHOLD;
	}


	
	public int getDopplerverschiebung(){
		return dopplerShift;
	}
	
	public int getCodeVerschiebung(){
		return codeShift;
	}

}
