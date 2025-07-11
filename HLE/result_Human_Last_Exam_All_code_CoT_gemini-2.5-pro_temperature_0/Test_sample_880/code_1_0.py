import pywt
import numpy as np

def solve_wavelet_choice():
    """
    Determines and explains the best mother wavelet for daily rainfall analysis.
    """
    # The chosen mother wavelet
    chosen_wavelet_name = 'db1'
    
    # Create the wavelet object
    try:
        wavelet = pywt.Wavelet(chosen_wavelet_name)
    except ValueError:
        print(f"Could not find the wavelet '{chosen_wavelet_name}'. Please ensure pywavelets is installed.")
        return

    print("### Analysis of the Best Mother Wavelet for Daily Rainfall Data ###\n")
    print("1. Signal Characteristics: Daily rainfall time series are characterized by abrupt changes and step-like behavior (i.e., sudden starts and stops of rain events), not smooth curves.")
    print("2. Analysis Goal: The goal is to preserve 'local details,' which means accurately capturing the timing of these sharp events.")
    print("3. Wavelet Choice Justification:")
    print(f"The '{wavelet.name}' (Daubechies1), also known as the Haar wavelet, is the best fit.")
    print("Its shape is a discontinuous step function, which naturally models the sharp jumps in rainfall data.")
    print("It has the most compact support, providing the best time localization to pinpoint exactly when rain events occur.\n")
    
    print(f"### Properties of the Chosen Wavelet: {wavelet.name} ###")
    print(f"Family Name: {wavelet.family_name}")
    print(f"Short Name: {wavelet.short_name}")
    print(f"Orthogonal: {wavelet.orthogonal}")
    print(f"Symmetric: {wavelet.symmetric}")
    
    # The filter coefficients can be considered the 'numbers in the final equation' that define the wavelet transform.
    print("\nDefining Equation (Decomposition Filter Coefficients):")
    
    # The decomposition low-pass filter (scaling function)
    dec_lo = wavelet.dec_lo
    print("Scaling function (low-pass filter) coefficients:")
    for num in dec_lo:
        print(f"{num:.8f}")

    # The decomposition high-pass filter (wavelet function)
    dec_hi = wavelet.dec_hi
    print("\nWavelet function (high-pass filter) coefficients:")
    for num in dec_hi:
        print(f"{num:.8f}")
        
    print("\nNote: The coefficients are sqrt(2)/2 for the low-pass filter and [-sqrt(2)/2, sqrt(2)/2] for the high-pass filter (in reverse order).")

solve_wavelet_choice()
print("\n<<<A>>>")