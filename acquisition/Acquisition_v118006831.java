
package gps.acquisition;

import cgra.pe.PETrigonometry;

public class Acquisition_v118006831 {

	// Static input parameters
	private static final int FREQ_STEP_COUNT = 11;
	private static final float FREQ_MOD_MIN = -5000;
	private static final float FREQ_MOD_MAX = 5000;
	private static final float FREQ_MOD_STEP = 1000;
	private static final float FREQ_SAMPLING = 400000;
	private static final float ACQ_THRESHOLD = 0.015f;
	private static final float FPI = (float) 3.14159265359;

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

	// Storage for DFT(code) results
	private final float[] code_r_dft;
	private final float[] code_i_dft;

	public Acquisition_v118006831(int nrOfSamples) {
		this.sCount = nrOfSamples;

		this.sample_r = new float[this.sCount];
		this.sample_i = new float[this.sCount];
		this.res_r_1 = new float[this.sCount];
		this.res_i_1 = new float[this.sCount];
		this.res_r_2 = new float[this.sCount];
		this.res_i_2 = new float[this.sCount];

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
		final float lFPI = FPI;
		final int sCount = this.sCount;
		float sumreal = 0;
		float sumimag = 0;
		final float constant = 2*lFPI/sCount;
		final float constant2 = 2 * lFPI/FREQ_SAMPLING;
		float pin = 0;
		
			
		// --- Prepare the code samples by transforming and complex conjugating them
		for (int k = 0; k < sCount; k++) { // For each output element
			pin = pin + sample_r[k] * sample_r[k] + sample_i[k] * sample_i[k];
			sumreal = 0;
			sumimag = 0;
			for (int t = 0; t < sCount; t++) { // For each input element
		//		float angle = 2 * lFPI * t * k / sCount;
				sumreal +=  code_r[t]  * PETrigonometry.cos( constant * t * k) + code_i[t]  * PETrigonometry.sin( constant * t * k);
				sumimag += -code_r[t]  * PETrigonometry.sin( constant * t * k) + code_i[t]  * PETrigonometry.cos( constant * t * k);
			}
			code_r_dft[k] 	=  sumreal;
			code_i_dft[k] 	= -sumimag;
		}
		

		// Compute the sin and cos, outside of the big loop. Define them as constants		
		final float xfd_cos2[][] = new float[FREQ_STEP_COUNT][sCount];
		final float xfd_sin2[][] = new float[FREQ_STEP_COUNT][sCount];	
		for (int j = 0; j < FREQ_STEP_COUNT; ++j) {
			final float fd = FREQ_MOD_MIN + FREQ_MOD_STEP * j;
			for (int i = 0; i < sCount; ++i)
			{
				xfd_cos2[j][i] =  PETrigonometry.cos(constant2 * fd * i );
				xfd_sin2[j][i] =  PETrigonometry.sin(constant2 * fd * i ) ;
			}
		}
		
		// --- Calculate Smax on the fly while performing all the other operations
		// OK
		float smax = 0;
		float fdAtMax = 0;
		int tauAtMax = 0;
		
		// --- Loop to visit all required frequency shifts (fd)		
		for (int j = 0; j < FREQ_STEP_COUNT; ++j)
		{
			// --- Value of fd for this loop iteration
			final float fd = FREQ_MOD_MIN + FREQ_MOD_STEP * j;

			// --- Wipe off the carrier (from the given samples) using the current fd
			for (int i = 0; i < sCount; ++i)
			{
				res_r_1[i] =  sample_r[i] * xfd_cos2[j][i] + sample_i[i] * xfd_sin2[j][i];
				res_i_1[i] = -sample_r[i] * xfd_sin2[j][i] + sample_i[i] * xfd_cos2[j][i];
			}

			// ... Samples in res1

			// --- Apply the DFT to the given samples which had the carrier wiped off
			// DFT FUNCTION INLINED
			for (int k = 0; k < sCount; k++) { // For each output element
				sumreal = 0;
				sumimag = 0;
				for (int t = 0; t < sCount; t++) { // For each input element
					// final float angle = constant * t * k;
					sumreal  +=  res_r_1[t] * PETrigonometry.cos(constant * t * k) + res_i_1[t] * PETrigonometry.sin(constant * t * k);
					sumimag  += -res_r_1[t] * PETrigonometry.sin(constant * t * k) + res_i_1[t] * PETrigonometry.cos(constant * t * k);
				}
				res_r_2[k] = sumreal;
				res_i_2[k] = sumimag;
			}

			// ... Samples in res2

			// --- Multiply both (samples and code) DFT results
			for (int i = 0; i < sCount; ++i) {
				res_r_1[i] = res_r_2[i] * code_r_dft[i] - res_i_2[i] * code_i_dft[i];
				res_i_1[i] = res_i_2[i] * code_r_dft[i] + res_r_2[i] * code_i_dft[i];
			}

			// ... Samples in res1
			// --- Apply the IDFT to retrieve the results
			// IDFT INLINED, combined with result search (smax, fdAtMax, tauAtMax)
			for (int k = 0; k < sCount; k++) { // For each output element
				sumreal = 0;
				sumimag = 0;
				
				for (int t = 0; t < sCount; t++) { // For each input element
					// final float angle = constant * t * k;
					sumreal += res_r_1[t] * PETrigonometry.cos(constant * t * k) - res_i_1[t] * PETrigonometry.sin(constant * t * k);
					sumimag += res_r_1[t] * PETrigonometry.sin(constant * t * k) + res_i_1[t] * PETrigonometry.cos(constant * t * k);
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
		// OK
		dopplerShift =  (int) fdAtMax;
		codeShift = tauAtMax;
		
		pin /= sCount;

		// --- Calculate the threshold by using the results
		// OK
		final float resultingThreshold = smax / pin;

		return resultingThreshold > ACQ_THRESHOLD;
	}


	
	public int getDopplerverschiebung(){
		return dopplerShift;
	}
	
	public int getCodeVerschiebung(){
		return codeShift;
	}

}