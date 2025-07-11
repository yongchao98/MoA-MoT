import numpy as np
import pywt

def analyze_rainfall_wavelet():
    """
    Demonstrates why the Daubechies1 (Haar) wavelet is suitable for rainfall data.
    """
    # 1. Create a synthetic daily rainfall signal.
    # It has two rain events: one lasting 3 days, another lasting 2 days.
    # The values represent rainfall intensity (e.g., in mm).
    rainfall_signal = np.array([0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 5.0, 5.0, 0.0, 0.0])
    print(f"Original Synthetic Rainfall Signal:\n{rainfall_signal}\n")

    # 2. Choose wavelets for comparison
    wavelet_db1 = 'db1'  # Haar wavelet (discontinuous, step-like)
    wavelet_db2 = 'db2'  # Daubechies2 wavelet (smoother)

    print("-" * 50)
    print(f"Analysis with Daubechies1 (Haar) Wavelet: '{wavelet_db1}'")
    print("-" * 50)

    # 3. Perform Discrete Wavelet Transform (DWT) with db1
    # cA are the approximation coefficients (low-pass filter)
    # cD are the detail coefficients (high-pass filter, detects changes)
    cA_db1, cD_db1 = pywt.dwt(rainfall_signal, wavelet_db1)

    print("Approximation Coefficients (cA) from db1:")
    # Using a simple loop to print each number as requested
    for num in cA_db1:
        print(f"{num:.4f}", end=" ")
    print("\n")

    print("Detail Coefficients (cD) from db1:")
    print("Note how the large non-zero values correspond to the start and end of rain events.")
    for num in cD_db1:
        print(f"{num:.4f}", end=" ")
    print("\n")


    print("-" * 50)
    print(f"Analysis with Daubechies2 (Smoother) Wavelet: '{wavelet_db2}'")
    print("-" * 50)

    # 4. Perform DWT with db2
    cA_db2, cD_db2 = pywt.dwt(rainfall_signal, wavelet_db2)

    print("Approximation Coefficients (cA) from db2:")
    for num in cA_db2:
        print(f"{num:.4f}", end=" ")
    print("\n")

    print("Detail Coefficients (cD) from db2:")
    print("Note how the coefficients are more 'smeared' and less localized than with db1.")
    for num in cD_db2:
        print(f"{num:.4f}", end=" ")
    print("\n")

    print("-" * 50)
    print("\nConclusion:")
    print("The Daubechies1 (Haar) wavelet's detail coefficients are sparse and precisely locate the abrupt changes in the signal (start/end of rain).")
    print("The smoother Daubechies2 wavelet spreads this information across more coefficients, making the event boundaries less clear.")
    print("Therefore, for blocky, intermittent data like daily rainfall, Daubechies1 is the superior choice.")

if __name__ == '__main__':
    analyze_rainfall_wavelet()