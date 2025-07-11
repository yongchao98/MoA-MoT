import numpy as np
import pywt

def demonstrate_wavelet_choice():
    """
    This function demonstrates why the Daubechies1 (Haar) wavelet is a good
    choice for rainfall-like signals with sharp transitions.
    """
    # 1. Create a synthetic signal representing a simple rainfall event.
    #    This signal is zero everywhere except for a sharp, blocky pulse.
    signal = np.zeros(32)
    signal[10:16] = 10.0  # A 6-day rainfall event with intensity 10

    print("--- Synthetic Rainfall Signal ---")
    print("A simple signal with a sharp start and end, similar to a rainfall event.")
    print(np.round(signal, 2))
    print("-" * 35)

    # 2. Decompose the signal using Daubechies1 (Haar) wavelet.
    #    The Haar wavelet is a step-function, ideal for capturing sharp jumps.
    wavelet_db1 = 'db1'
    coeffs_db1 = pywt.dwt(signal, wavelet_db1, mode='symmetric')
    cA_db1, cD_db1 = coeffs_db1

    print(f"--- Analysis with {wavelet_db1} (Haar wavelet) ---")
    print("The detail coefficients (cD) show a very localized response.")
    print("The sharp start/end of the rain event are captured by two strong coefficients.")
    # The detail coefficients capture the change/difference in the signal.
    # The Haar wavelet essentially calculates differences between adjacent pairs of averages.
    # The large values correspond to where the signal jumps from 0 to 10 and back.
    print(np.round(cD_db1, 2))
    print("-" * 35)

    # 3. Decompose the signal using a smoother wavelet, Daubechies2.
    #    This will show how smoother wavelets can "smear" the event.
    wavelet_db2 = 'db2'
    coeffs_db2 = pywt.dwt(signal, wavelet_db2, mode='symmetric')
    cA_db2, cD_db2 = coeffs_db2

    print(f"--- Analysis with {wavelet_db2} (smoother wavelet) ---")
    print("The detail coefficients (cD) are less localized (smeared).")
    print("The energy of the sharp jump is spread across multiple coefficients.")
    print(np.round(cD_db2, 2))
    print("-" * 35)

    print("\nConclusion:")
    print("The Daubechies1 wavelet provides a more compact and interpretable representation")
    print("of sharp, intermittent events like daily rainfall, making it the best fit.")

# Run the demonstration
demonstrate_wavelet_choice()
<<<A>>>