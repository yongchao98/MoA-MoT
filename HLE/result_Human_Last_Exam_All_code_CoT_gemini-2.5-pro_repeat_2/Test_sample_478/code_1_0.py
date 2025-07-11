import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio at a specific location.
    """
    # Step 1: Define given values and assumptions.
    # L_signal is the sound level of the people at their location.
    L_signal = 75.0  # dB

    # We assume the 75 dB level is the source strength at a reference distance of 1m.
    L_people_source = 75.0 # dB at 1m
    r_people_source = 1.0  # m

    # Distances for the final SNR calculation.
    r_train_final = 30.0  # m
    r_construction_final = 50.0  # m

    # Measurement A data: used to find the train's noise level.
    L_total_A = 100.0  # dB
    r_train_A = 10.0   # m
    r_people_A = 20.0  # m

    # Measurement B data: used to find the construction's noise level.
    L_total_B = 115.0  # dB
    r_construction_B = 20.0 # m
    r_people_B = 30.0   # m

    # Step 2: Determine the sound level of the train at 10m (L_train_A).
    # First, find the sound level of the people at measurement point A.
    L_people_at_A = L_people_source - 20 * math.log10(r_people_A / r_people_source)
    
    # Convert levels to relative intensities to subtract the people's contribution.
    # Note: I_rel = 10^(L/10)
    I_total_A_rel = 10**(L_total_A / 10)
    I_people_A_rel = 10**(L_people_at_A / 10)
    
    # The train's intensity is the total minus the people's.
    I_train_A_rel = I_total_A_rel - I_people_A_rel
    
    # Convert the train's intensity back to decibels.
    L_train_A = 10 * math.log10(I_train_A_rel)

    # Step 3: Determine the sound level of the construction at 20m (L_construction_B).
    # First, find the sound level of the people at measurement point B.
    L_people_at_B = L_people_source - 20 * math.log10(r_people_B / r_people_source)

    # Convert levels to relative intensities.
    I_total_B_rel = 10**(L_total_B / 10)
    I_people_B_rel = 10**(L_people_at_B / 10)

    # The construction's intensity is the total minus the people's.
    I_construction_B_rel = I_total_B_rel - I_people_B_rel

    # Convert the construction's intensity back to decibels.
    L_construction_B = 10 * math.log10(I_construction_B_rel)

    # Step 4: Calculate the noise levels at the final location (where the people are).
    # Use the inverse square law for decibels: L2 = L1 - 20 * log10(r2 / r1)
    L_train_final = L_train_A - 20 * math.log10(r_train_final / r_train_A)
    L_construction_final = L_construction_B - 20 * math.log10(r_construction_final / r_construction_B)

    # Step 5: Combine the two noise sources to get the total noise level.
    # Convert individual noise levels to intensities, add them, then convert back to dB.
    I_train_final_rel = 10**(L_train_final / 10)
    I_construction_final_rel = 10**(L_construction_final / 10)
    I_noise_total_rel = I_train_final_rel + I_construction_final_rel
    L_noise_total = 10 * math.log10(I_noise_total_rel)

    # Step 6: Calculate the Signal-to-Noise Ratio (SNR).
    snr = L_signal - L_noise_total

    # Print the final equation with all the calculated numbers.
    print(f"Signal Level: {L_signal:.2f} dB")
    print(f"Noise from Train at {r_train_final}m: {L_train_final:.2f} dB")
    print(f"Noise from Construction at {r_construction_final}m: {L_construction_final:.2f} dB")
    print(f"Total Noise Level: {L_noise_total:.2f} dB")
    print("\n--- Final Calculation ---")
    print(f"SNR = Signal Level - Total Noise Level")
    print(f"SNR = {L_signal:.2f} dB - {L_noise_total:.2f} dB = {snr:.2f} dB")

solve_snr()
<<<D>>>