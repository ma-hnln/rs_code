package gps.acquisition;


public class Acquisition {

	// Static input parameters
	private static final int FREQ_STEP_COUNT = 11;
	private static final float FREQ_MOD_MIN = -5000;
	private static final float FREQ_MOD_MAX = 5000;
	private static final float FREQ_MOD_STEP = 1000;
	private static final float FREQ_SAMPLING = 400000;
	private static final float ACQ_THRESHOLD = 0.015f;

	// Results
	private int dopplerShift = 0;
	private int codeShift = 0;

	// Dynamic input parameters
	private final int sCount;

	// Sample storage, currently only for unmodified samples
	private final float[] sample_r;
	private final float[] sample_i;
	private int sampleWriteIndex = 0;

	// Code storage
	private final float[] code_r;
	private final float[] code_i;
	private int codeWriteIndex = 0;

	public Acquisition(int nrOfSamples) {
		this.sCount = nrOfSamples;
		this.sample_r = new float[this.sCount];
		this.sample_i = new float[this.sCount];
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
		++sampleWriteIndex;
	}
	
	public boolean startAcquisition(){
		// Code for one (un-shifted) sample pack, see sample_r and sample_i
		final float FPI = (float) Math.PI;

		// Wipe off the carrier (from the given samples)
		// TODO: for all fd
		{
			float fd = 0; // We currently use the un-shifted samples

			for (int i = 0; i < sCount; ++i)
			{
				final float angle = 2 * FPI * fd * i / FREQ_SAMPLING;
				final float sined = (float) Math.sin(angle);
				final float cosed = (float) Math.cos(angle);

				final float s_r =  sample_r[i] * cosed + sample_i[i] * sined;
				final float s_i = -sample_r[i] * sined + sample_i[i] * cosed;

				sample_r[i] = s_r;
				sample_i[i] = s_i;
			}
		}

		// Apply DFT to the given code samples and complex conjugate the result
		// -> What about the transpose operation from the slides?
		// TODO: for all fd
		float[] code_r_dft = new float[sCount];
		float[] code_i_dft = new float[sCount];

		dft(code_r, code_i, code_r_dft, code_i_dft);

		for (int i = 0; i < sCount; ++i)
			code_i_dft[i] *= -1;

		// Apply DFT to the given samples which had the carrier wiped off
		// TODO: for all fd
		float[] sample_r_dft = new float[sCount];
		float[] sample_i_dft = new float[sCount];

		dft(sample_r, sample_i, sample_r_dft, sample_i_dft);

		// Multiply both DFT results
		// TODO: how is this multiplication performed? also, for all fd
		float[] sample_code_mul_r = new float[sCount];
		float[] sample_code_mul_i = new float[sCount];

		// Apply IDFT (normalise the result at the same time?)
		// TODO: for all fd

		float[] res_idft_r = new float[sCount];
		float[] res_idft_i = new float[sCount];

		idft(sample_code_mul_r, sample_code_mul_i, res_idft_r, res_idft_i);

		// Calculate Smax of the resulting vector (will be the matrix later on)
		float smax = 0;
		int fd = 0;
		int tau = 0;

		// TODO: for all fd
		for (int i = 0; i < sCount; ++i) {
			final float abs_sqr = res_idft_r[i] * res_idft_r[i] + res_idft_i[i] * res_idft_i[i];
			if (smax < abs_sqr) {
				smax = abs_sqr;
				fd = 0; // To be replaced by the actual loop variant
				tau = i;
			}
		}

		dopplerShift = fd;
		codeShift = tau;

		// Calculate the signals' power (do not forget the normalisation)
		float pin = 0;
		for (int i = 0; i < sCount; ++i) {
			final float abs_sqr = sample_r[i] * sample_r[i] + sample_i[i] * sample_i[i];
			pin += abs_sqr;
		}
		pin /= sCount;

		// Calculate the threshold from these values
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
		final int n = in_r.length;
		for (int k = 0; k < n; k++) { // For each output element
			float sumreal = 0;
			float sumimag = 0;
			for (int t = 0; t < n; t++) { // For each input element
				double angle = 2 * Math.PI * t * k / n;
				sumreal +=  in_r[t] * Math.cos(angle) + in_i[t] * Math.sin(angle);
				sumimag += -in_r[t] * Math.sin(angle) + in_i[t] * Math.cos(angle);
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
		final int n = in_r.length;
		for (int k = 0; k < n; k++) { // For each output element
			float sumreal = 0;
			float sumimag = 0;
			for (int t = 0; t < n; t++) { // For each input element
				double angle = 2 * Math.PI * t * k / n;
				sumreal +=  in_r[t] * Math.sin(angle) + in_i[t] * Math.cos(angle);
				sumimag += -in_r[t] * Math.cos(angle) + in_i[t] * Math.sin(angle);
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
