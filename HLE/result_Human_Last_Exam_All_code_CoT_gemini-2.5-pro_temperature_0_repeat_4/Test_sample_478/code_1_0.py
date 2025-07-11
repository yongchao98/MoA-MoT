import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio at a specific location based on given sound sources and levels.
    """
    # Step 1: Define the given values from the problem description.
    # Signal level at the location of interest.
    L_signal = 75.0  # dB

    # Known noise level of the train at a reference distance.
    L_train_ref = 100.0  # dB
    dist_train_ref = 10.0  # meters

    # Known noise level of the construction at a reference distance.
    L_constr_ref = 115.0  # dB
    dist_constr_ref = 20.0  # meters

    # Distances from the noise sources to the location of the people.
    dist_train_to_people = 30.0  # meters
    dist_constr_to_people = 50.0  # meters

    # Step 2: Calculate the sound level of each noise source at the people's location.
    # We use the formula: L2 = L1 - 20 * log10(r2 / r1)

    # For the train:
    L_train_at_people = L_train_ref - 20 * math.log10(dist_train_to_people / dist_train_ref)

    # For the construction:
    L_constr_at_people = L_constr_ref - 20 * math.log10(dist_constr_to_people / dist_constr_ref)

    # Step 3: Combine the two noise levels.
    # This requires converting dB to intensities, adding them, and converting back to dB.
    # I_relative = 10^(L/10)
    I_train_relative = 10**(L_train_at_people / 10)
    I_constr_relative = 10**(L_constr_at_people / 10)

    # Total relative intensity of the noise
    I_noise_total_relative = I_train_relative + I_constr_relative

    # Convert total relative intensity back to a dB level
    # L = 10 * log10(I_relative)
    L_noise_total = 10 * math.log10(I_noise_total_relative)

    # Step 4: Calculate the final Signal-to-Noise Ratio (SNR).
    # SNR (dB) = Signal Level (dB) - Total Noise Level (dB)
    snr = L_signal - L_noise_total

    # Print the final equation with all the calculated numbers.
    print("The Signal-to-Noise Ratio (SNR) is calculated as: Signal Level - Total Noise Level.")
    print("\nFirst, the total noise level is found by combining the noise from the train and construction:")
    print(f"Noise Level = 10 * log10( 10^(Train Level/10) + 10^(Construction Level/10) )")
    print(f"Noise Level = 10 * log10( 10^({L_train_at_people:.2f}/10) + 10^({L_constr_at_people:.2f}/10) ) = {L_noise_total:.2f} dB")
    
    print("\nNow, we can calculate the SNR:")
    print(f"SNR = Signal Level - Total Noise Level")
    print(f"SNR = {L_signal:.2f} dB - {L_noise_total:.2f} dB = {snr:.2f} dB")

solve_snr()
<<<D>>>