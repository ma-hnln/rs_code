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
	private int sCount = 0;

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
		// TODO
		{
			float fd = 0; // We use the un-shifted samples

			for (int i = 0; i < sCount; ++i) {
				final float angle = 2 * FPI * fd * i / FREQ_SAMPLING;
				final float s_r = sample_r[i] * (float) Math.cos(angle) + sample_i[i] * (float) Math.sin(angle);
			}
		}

		// Apply DFT to the given code samples and complex conjugate the result
		// -> What about the transpose operation from the slides?
		// TODO

		// Apply DFT to the given samples which had the carrier wiped off
		// TODO

		// Multiply both DFT results
		// TODO

		// Apply IDFT (normalise the result at the same time?)
		// TODO

		// Calculate Smax of the resulting vector (will be the matrix later on)
		float smax = 0;
		int fd = 0;
		int tau = 0;

		// TODO

		dopplerShift = fd;
		codeShift = tau;

		// Calculate the signals' power (do not forget the normalisation)
		final float pin = 0;

		// TODO

		// Calculate the threshold from these values
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
