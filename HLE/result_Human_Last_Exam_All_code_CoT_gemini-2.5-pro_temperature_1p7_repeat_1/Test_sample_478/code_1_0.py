import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio at the location of the people.
    """
    # Signal level at the location of the people (given)
    L_signal = 75.0  # in dB

    # Measurement 1 data
    L_total_M1 = 100.0  # in dB
    r_train_M1 = 10.0   # in meters

    # Measurement 2 data
    L_total_M2 = 115.0  # in dB
    r_construction_M2 = 20.0  # in meters
    
    # Distances from noise sources to the people's location
    r_train_to_people = 30.0  # in meters
    r_construction_to_people = 50.0 # in meters

    # Step 1: Calculate the source strength (at 1m) of the train.
    # We assume the 100 dB at 10m is dominated by the train.
    # L_1m = L_10m - 20 * log10(1/10) = L_10m + 20 * log10(10)
    L_train_1m = L_total_M1 + 20 * math.log10(r_train_M1)

    # Step 2: Calculate the source strength (at 1m) of the construction.
    # We assume the 115 dB at 20m is dominated by the construction.
    # L_1m = L_20m - 20 * log10(1/20) = L_20m + 20 * log10(20)
    L_construction_1m = L_total_M2 + 20 * math.log10(r_construction_M2)

    # Step 3: Calculate the noise level from the train at the people's location (30m).
    # L_30m = L_1m - 20 * log10(30/1)
    L_train_at_people = L_train_1m - 20 * math.log10(r_train_to_people)

    # Step 4: Calculate the noise level from the construction at the people's location (50m).
    # L_50m = L_1m - 20 * log10(50/1)
    L_construction_at_people = L_construction_1m - 20 * math.log10(r_construction_to_people)
    
    # Step 5: Combine the two noise levels to get the total noise at the people's location.
    # L_total = 10 * log10(10^(L_A/10) + 10^(L_B/10))
    intensity_train = 10**(L_train_at_people / 10)
    intensity_construction = 10**(L_construction_at_people / 10)
    L_noise_total = 10 * math.log10(intensity_train + intensity_construction)

    # Step 6: Calculate the Signal-to-Noise Ratio (SNR).
    snr = L_signal - L_noise_total

    print("Step-by-step calculation:")
    print(f"Signal Level (S): {L_signal:.2f} dB")
    print(f"Noise level from train at people's location: {L_train_at_people:.2f} dB")
    print(f"Noise level from construction at people's location: {L_construction_at_people:.2f} dB")
    print(f"Total Noise Level (N): {L_noise_total:.2f} dB")
    print("\nFinal Equation (SNR = S - N):")
    # Outputting each number in the final equation
    print(f"SNR = {L_signal:.2f} dB - {L_noise_total:.2f} dB = {snr:.2f} dB")
    
solve_snr()
<<<D>>>