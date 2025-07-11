import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio at the location of the people.
    """
    # Step 1: Define the known values from the problem description.
    # We approximate the source levels by ignoring the minor contribution from the people's noise at the measurement points.
    l_train_at_10m = 100  # dB
    l_construction_at_20m = 115  # dB
    
    # Distances of noise sources from the people's location
    dist_train_to_people = 30  # meters
    dist_construction_to_people = 50  # meters
    
    # Signal level at the people's location
    l_signal = 75  # dB

    # Step 2: Calculate the noise level from each source at the people's location.
    # Formula: L2 = L1 - 20 * log10(r2 / r1)
    
    # Noise from the train at the people's location (30m away)
    l_train_at_people = l_train_at_10m - 20 * math.log10(dist_train_to_people / 10)
    
    # Noise from the construction at the people's location (50m away)
    l_construction_at_people = l_construction_at_20m - 20 * math.log10(dist_construction_to_people / 20)

    # Step 3: Combine the two noise levels to get the total noise.
    # Convert dB to intensity, add, and convert back to dB.
    # L_total = 10 * log10(10^(L1/10) + 10^(L2/10))
    intensity_train = 10**(l_train_at_people / 10)
    intensity_construction = 10**(l_construction_at_people / 10)
    l_total_noise = 10 * math.log10(intensity_train + intensity_construction)

    # Step 4: Calculate the Signal-to-Noise Ratio (SNR).
    # SNR = Signal (dB) - Noise (dB)
    snr = l_signal - l_total_noise

    # Print the breakdown of the final calculation
    print(f"The final calculation is SNR = Signal - Total Noise")
    print(f"Signal Level = {l_signal:.2f} dB")
    print(f"Calculated noise from train at people's location = {l_train_at_people:.2f} dB")
    print(f"Calculated noise from construction at people's location = {l_construction_at_people:.2f} dB")
    print(f"Total Noise Level = 10 * log10(10^({l_train_at_people:.2f}/10) + 10^({l_construction_at_people:.2f}/10)) = {l_total_noise:.2f} dB")
    print(f"\nFinal Equation:")
    print(f"SNR = {l_signal:.2f} dB - {l_total_noise:.2f} dB")
    print(f"SNR = {snr:.2f} dB")

solve_snr()