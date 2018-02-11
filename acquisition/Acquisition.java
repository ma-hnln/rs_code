package gps.acquisition;

import cgra.pe.PETrigonometry;

public class Acquisition {

	// Static input parameters
	private static final int FREQ_STEP_COUNT = 11;
	private static final float FREQ_MOD_MIN = -5000;
	private static final float FREQ_MOD_MAX = 5000;
	private static final float FREQ_MOD_STEP = 1000;
	private static final float FREQ_SAMPLING = 400000;
	private static final float ACQ_THRESHOLD = 0.015f;
	private static final float FPI = (float) Math.PI;

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
	private final float[] sample_dft_r;
	private final float[] sample_dft_i;

	private final float[] res_r_1;
	private final float[] res_i_1;

	// Code storage (not transformed)
	private final float[] code_r;
	private final float[] code_i;
	private int codeWriteIndex = 0;

	// Storage for DFT(code) results
	private final float[] code_r_dft;
	private final float[] code_i_dft;

	public Acquisition(int nrOfSamples) {
		this.sCount = nrOfSamples;

		this.sample_r = new float[this.sCount];
		this.sample_i = new float[this.sCount];

		this.sample_dft_r = new float[this.sCount];
		this.sample_dft_i = new float[this.sCount];
		this.res_r_1 = new float[this.sCount];
		this.res_i_1 = new float[this.sCount];

		this.code_r = new float[this.sCount];
		this.code_i = new float[this.sCount];
		this.code_r_dft = new float[this.sCount];
		this.code_i_dft = new float[this.sCount];
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
		
		final float twoPIdivByN = 2 * FPI / sCount;
		float sumreal = 0;
		float sumimag = 0;
		
		// --- Prepare the code samples by transforming and complex conjugating them
		for (int k = 0; k < sCount; k++) { // For each output element
			sumreal = 0;
			sumimag = 0;
			for (int t = 0; t < sCount; t++) { // For each input element
				sumreal +=  code_r[t] * PETrigonometry.cos(twoPIdivByN * t * k) + code_i[t] * PETrigonometry.sin(twoPIdivByN * t * k);
				sumimag += -code_r[t] * PETrigonometry.sin(twoPIdivByN * t * k) + code_i[t] * PETrigonometry.cos(twoPIdivByN * t * k);
			}
			code_r_dft[k] =  sumreal;
			code_i_dft[k] = -sumimag;
		}

		// --- Prepare the samples by by transforming them
		// We do this with the raw samples, no carrier wipe off is performed at this point.
		// The wipe off is happening by applying it in the frequency domain.
		
		// --- Also calculate the signals' power (do not forget the normalisation)
		float pin = 0;
		for (int k = 0; k < sCount; k++) { // For each output element
			pin += sample_r[k] * sample_r[k] + sample_i[k] * sample_i[k];
			sumreal = 0;
			sumimag = 0;
			for (int t = 0; t < sCount; t++) { // For each input element
				sumreal +=  sample_r[t] * PETrigonometry.cos(twoPIdivByN * t * k) + sample_i[t] * PETrigonometry.sin(twoPIdivByN * t * k);
				sumimag += -sample_r[t] * PETrigonometry.sin(twoPIdivByN * t * k) + sample_i[t] * PETrigonometry.cos(twoPIdivByN * t * k);
			}
			sample_dft_r[k] = sumreal;
			sample_dft_i[k] = sumimag;
		}
		pin /= sCount;
		
		// --- Calculate Smax on the fly while performing all the other operations
		float smax = 0;
		float fdAtMax = 0;
		int tauAtMax = 0;

		// --- Loop to visit all required frequency shifts (fd)
		for (int j = 0; j < FREQ_STEP_COUNT; ++j)
		{
			// --- Value of fd for this loop iteration
			final float fd = FREQ_MOD_MIN + FREQ_MOD_STEP * j;

			// --- Wipe off the carrier (from the given samples) using the current fd
			// We do this in the frequency domain by applying the frequency shift
			// directly to the transformed samples: 
			//
			// x(n)*e^(j*2*PI*n*i/N) , with i = -N*fd/fs
			// ->
			// X((k-i) mod N) 
			//
			// We omit the '-' of i as we will just add it to k. 
			{
				// This is going to be a problem for lower sample counts...
				final int i = (sCount * (int) fd) / (int) FREQ_SAMPLING;
				// final int i =  Math.round(((float)sCount * fd) / FREQ_SAMPLING);
				final int read_start = (sCount + i) % sCount;
				final int first_write_end = sCount - read_start;
				
				int read_index = read_start;
				int k = 0;
				
				// --- Directly multiply both (samples and code) DFT results
				for (; k < first_write_end; ++k)
				{
					res_r_1[k] = sample_dft_r[read_index] * code_r_dft[k] - sample_dft_i[read_index] * code_i_dft[k]; 
					res_i_1[k] = sample_dft_i[read_index] * code_r_dft[k] + sample_dft_r[read_index] * code_i_dft[k];
					++read_index;
				}
				
				read_index = 0; // Wrapping around to zero
				for (; k < sCount; ++k)
				{
					res_r_1[k] = sample_dft_r[read_index] * code_r_dft[k] - sample_dft_i[read_index] * code_i_dft[k]; 
					res_i_1[k] = sample_dft_i[read_index] * code_r_dft[k] + sample_dft_r[read_index] * code_i_dft[k];
					++read_index;
				}
			}

			// ... Samples in res1

			// --- Apply the IDFT to retrieve the results
			for (int k = 0; k < sCount; k++) { // For each output element
				sumreal = 0;
				sumimag = 0;
				
				for (int t = 0; t < sCount; t++) { // For each input element
					// final float angle = constant * t * k;
					sumreal += res_r_1[t] * PETrigonometry.cos(twoPIdivByN * t * k) - res_i_1[t] * PETrigonometry.sin(twoPIdivByN * t * k);
					sumimag += res_r_1[t] * PETrigonometry.sin(twoPIdivByN * t * k) + res_i_1[t] * PETrigonometry.cos(twoPIdivByN * t * k);
				}
				
				final float res_r = sumreal / sCount;
				final float res_i = sumimag / sCount;
				
				// --- Search the maximum in the resulting vector
				if (smax < res_r * res_r + res_i * res_i) {
					smax = res_r * res_r + res_i * res_i;
					fdAtMax = fd;
					tauAtMax = k;
				}
			}
		}

		// --- Save results
		dopplerShift =  (int) fdAtMax;
		codeShift = tauAtMax;

		// --- Calculate the threshold by using the results
		final float resultingThreshold = smax / pin;

		return resultingThreshold > ACQ_THRESHOLD;
	}

	/**
	 * Perform the DFT for the given array of complex float values.
	 * Assumes that all arrays are of equal length N.
	 * Initial version of tho code taken from:
	 * https://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
	 *
	 * @param in_r real part of the input values
	 * @param in_i imaginary part of the input values
	 * @param out_r real part of the output values
	 * @param out_i imaginary part of the output values
	 */
	static private void dft(float[] in_r , float[] in_i,
					 		float[] out_r, float[] out_i)
	{
		final float lFPI = FPI;
		final int n = in_r.length;

		for (int k = 0; k < n; k++) { // For each output element
			float sumreal = 0;
			float sumimag = 0;
			for (int t = 0; t < n; t++) { // For each input element
				float angle = 2 * lFPI * t * k / n;
				sumreal +=  in_r[t] * PETrigonometry.cos(angle) + in_i[t] * PETrigonometry.sin(angle);
				sumimag += -in_r[t] * PETrigonometry.sin(angle) + in_i[t] * PETrigonometry.cos(angle);
			}
			out_r[k] = sumreal;
			out_i[k] = sumimag;
		}
	}

	/**
	 * Perform the IDFT for the given array of complex float values.
	 * Assumes that all arrays are of equal length N.
	 *
	 * @param in_r real part of the input values
	 * @param in_i imaginary part of the input values
	 * @param out_r real part of the output values
	 * @param out_i imaginary part of the output values
	 */
	static private void idft(float[] in_r , float[] in_i,
							 float[] out_r, float[] out_i)
	{
		final float lFPI = FPI;
		final int n = in_r.length;

		for (int k = 0; k < n; k++) { // For each output element
			float sumreal = 0;
			float sumimag = 0;
			for (int t = 0; t < n; t++) { // For each input element
				float angle = 2 * lFPI * t * k / n;
				sumreal += in_r[t] * PETrigonometry.cos(angle) - in_i[t] * PETrigonometry.sin(angle);
				sumimag += in_r[t] * PETrigonometry.sin(angle) + in_i[t] * PETrigonometry.cos(angle);
			}
			out_r[k] = sumreal / n;
			out_i[k] = sumimag / n;
		}
	}
	
	public int getDopplerverschiebung(){
		return dopplerShift;
	}
	
	public int getCodeVerschiebung(){
		return codeShift;
	}

}
