import math

def calculate_snr():
    """
    This script calculates the signal-to-noise ratio (SNR) based on the problem description.
    """
    # --- Step 1: Define known values from the problem ---
    S = 75.0  # Signal level of the people in dB

    # Train noise data
    l_train_r1 = 100.0  # dB
    r_train_1 = 30.0    # meters (distance to people)
    r_train_ref = 10.0  # meters (reference distance)

    # Construction noise data
    l_const_r1 = 115.0  # dB
    r_const_1 = 50.0    # meters (distance to people)
    r_const_ref = 20.0  # meters (reference distance)
    
    # --- Step 2: Calculate the noise level from the train at the people's location ---
    # Using the formula: L2 = L1 - 20 * log10(r2 / r1)
    l_train_at_people = l_train_r1 - 20 * math.log10(r_train_1 / r_train_ref)
    
    # --- Step 3: Calculate the noise level from construction at the people's location ---
    l_const_at_people = l_const_r1 - 20 * math.log10(r_const_1 / r_const_ref)

    # --- Step 4: Calculate the total noise level (N) by combining the two sources ---
    # Convert dB to intensity, sum them, and convert back to dB.
    # I_total = 10^(L_train/10) + 10^(L_const/10)
    # N = 10 * log10(I_total)
    intensity_train = 10**(l_train_at_people / 10)
    intensity_const = 10**(l_const_at_people / 10)
    total_intensity = intensity_train + intensity_const
    N = 10 * math.log10(total_intensity)

    # --- Step 5: Calculate the Signal-to-Noise Ratio (SNR) ---
    snr = S - N

    # --- Step 6: Print the results ---
    print(f"The signal level (S) is given as: {S:.2f} dB")
    print("-" * 30)
    print(f"Calculating noise at the people's location:")
    print(f"Train noise at {r_train_1}m = {l_train_at_people:.2f} dB")
    print(f"Construction noise at {r_const_1}m = {l_const_at_people:.2f} dB")
    print(f"The total combined noise level (N) is: {N:.2f} dB")
    print("-" * 30)
    print("The final SNR is calculated as S - N:")
    print(f"{S:.2f} dB - {N:.2f} dB = {snr:.2f} dB")

# Execute the calculation
calculate_snr()
<<<D>>>