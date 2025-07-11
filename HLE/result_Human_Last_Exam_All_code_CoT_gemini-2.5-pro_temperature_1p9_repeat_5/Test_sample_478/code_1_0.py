import math

def solve_snr():
    """
    This function calculates the signal-to-noise ratio (SNR) at the specified location.
    """
    # Step 1: Define signal level and initial measurement parameters.
    # The signal (S) is the sound of the people talking.
    signal_db = 75.0

    # Approximated noise levels from the problem description.
    # The contribution of the people's talking is negligible at the measurement points.
    l_train_at_10m = 100.0
    l_construction_at_20m = 115.0
    
    # Distances for calculating noise at the people's location.
    dist_train_to_people = 30.0
    dist_construction_to_people = 50.0
    
    # Reference distances from the measurements.
    ref_dist_train = 10.0
    ref_dist_construction = 20.0

    # Step 2: Calculate noise level from each source at the people's location.
    # Sound attenuation formula: L2 = L1 - 20 * log10(r2 / r1)

    # Noise from the train at the people's location
    l_train_at_people = l_train_at_10m - 20 * math.log10(dist_train_to_people / ref_dist_train)

    # Noise from the construction at the people's location
    l_construction_at_people = l_construction_at_20m - 20 * math.log10(dist_construction_to_people / ref_dist_construction)
    
    # Step 3: Calculate the total combined noise level (N).
    # Convert dB to intensity, add intensities, and convert back to dB.
    # I = 10^(L/10)
    intensity_train = 10**(l_train_at_people / 10)
    intensity_construction = 10**(l_construction_at_people / 10)
    
    total_intensity_noise = intensity_train + intensity_construction
    
    # L_total = 10 * log10(I_total)
    total_noise_db = 10 * math.log10(total_intensity_noise)

    # Step 4: Calculate the final Signal-to-Noise Ratio (SNR).
    # SNR = S - N
    snr = signal_db - total_noise_db
    
    # Print the final equation as requested
    print(f"Signal Level (S): {signal_db:.2f} dB")
    print(f"Noise from Train: {l_train_at_people:.2f} dB")
    print(f"Noise from Construction: {l_construction_at_people:.2f} dB")
    print(f"Total Noise Level (N): {total_noise_db:.2f} dB")
    print(f"Final Equation: SNR = {signal_db:.2f} dB - {total_noise_db:.2f} dB")
    print(f"Signal-to-Noise Ratio (SNR): {snr:.2f} dB")

solve_snr()
<<<D>>>