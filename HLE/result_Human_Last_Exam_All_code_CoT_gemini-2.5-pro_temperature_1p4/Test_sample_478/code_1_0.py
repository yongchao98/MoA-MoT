import math

def calculate_snr():
    """
    Calculates the signal-to-noise ratio at a specific location based on given sound levels and distances.
    """

    # Step 1: Define known values from the problem statement.
    L_signal = 75.0  # dB, Signal level at the location of the people.

    # Information about the noise sources based on measurements.
    # We approximate the source levels by ignoring the much quieter sound from the people at the measurement points.
    L_train_at_10m = 100.0  # dB
    L_const_at_20m = 115.0  # dB

    # Distances from the noise sources to the location of the people (where SNR is calculated).
    r_train_final = 30.0  # meters
    r_const_final = 50.0  # meters
    
    # Reference distances for the noise source measurements.
    r_train_ref = 10.0 # meters
    r_const_ref = 20.0 # meters

    # Step 2: Calculate the noise level from each source at the people's location.
    # Formula: L2 = L1 - 20 * log10(r2 / r1)
    L_train_final = L_train_at_10m - 20 * math.log10(r_train_final / r_train_ref)
    L_const_final = L_const_at_20m - 20 * math.log10(r_const_final / r_const_ref)

    # Step 3: Combine the noise levels.
    # Convert dB to intensity ratio, add them, then convert back to dB.
    # Intensity_Ratio = 10^(L/10)
    intensity_ratio_train = 10**(L_train_final / 10)
    intensity_ratio_const = 10**(L_const_final / 10)

    # Total noise intensity is the sum of individual intensities.
    total_intensity_ratio = intensity_ratio_train + intensity_ratio_const

    # Convert total intensity back to dB: L = 10 * log10(Intensity_Ratio)
    L_total_noise = 10 * math.log10(total_intensity_ratio)

    # Step 4: Calculate the Signal-to-Noise Ratio (SNR).
    snr = L_signal - L_total_noise
    
    # Output the final equation and the result.
    print(f"The signal level is {L_signal:.2f} dB.")
    print(f"The noise level from the train at the location is {L_train_final:.2f} dB.")
    print(f"The noise level from the construction at the location is {L_const_final:.2f} dB.")
    print(f"The total combined noise level is {L_total_noise:.2f} dB.")
    print("\nThe final SNR is calculated as Signal Level - Total Noise Level.")
    print(f"Final Equation: SNR = {L_signal:.2f} dB - {L_total_noise:.2f} dB")
    print(f"\nCalculated Signal-to-Noise Ratio (SNR): {snr:.2f} dB")

# Run the calculation.
calculate_snr()