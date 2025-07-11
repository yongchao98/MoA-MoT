import math

def solve_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define the signal level at the people's location.
    L_S = 75  # dB

    # Step 2: Use the provided measurements to establish the strength of the noise sources.
    # As the train and construction noise levels are much higher than the people's talking,
    # we can approximate that the measured levels are dominated by the main noise sources.
    # Noise level of the train at 10 meters.
    L_train_10m = 100  # dB
    # Noise level of the construction at 20 meters.
    L_construction_20m = 115  # dB

    # Step 3: Calculate the noise level from each source at the people's location
    # using the formula: L2 = L1 - 20 * log10(r2 / r1)

    # People are 30m from the train
    r_train_initial = 10
    r_train_final = 30
    L_train_at_people = L_train_10m - 20 * math.log10(r_train_final / r_train_initial)

    # People are 50m from the construction
    r_construction_initial = 20
    r_construction_final = 50
    L_construction_at_people = L_construction_20m - 20 * math.log10(r_construction_final / r_construction_initial)

    # Step 4: Calculate the total noise level by adding the intensities.
    # Convert dB to intensity: I_norm = 10^(L/10)
    intensity_train = 10**(L_train_at_people / 10)
    intensity_construction = 10**(L_construction_at_people / 10)
    
    # Add intensities
    total_noise_intensity = intensity_train + intensity_construction
    
    # Convert total intensity back to dB: L = 10 * log10(I_norm)
    L_N_total = 10 * math.log10(total_noise_intensity)

    # Step 5: Calculate the Signal-to-Noise Ratio (SNR).
    snr = L_S - L_N_total

    # Step 6: Print the components of the final calculation and the result.
    print(f"The signal level (L_S) is {L_S} dB.")
    print(f"The total combined noise level (L_N) at the location is {L_N_total:.2f} dB.")
    print("\nThe final equation for Signal-to-Noise Ratio (SNR) is: L_S - L_N")
    print(f"SNR = {L_S:.2f} dB - {L_N_total:.2f} dB = {snr:.2f} dB")

solve_snr()
<<<D>>>