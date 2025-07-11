import math

def calculate_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # --- Step 1: Define constants and assumptions from the problem ---
    # Signal level at the location of the people
    L_signal = 75  # dB

    # Assumption: The 75 dB level of the people is measured at a reference distance of 1m.
    L_people_ref = 75  # dB
    d_people_ref = 1    # m

    # Distances for the two measurement points
    d_train_m1 = 10
    d_people_m1 = 20
    L_total_m1 = 100

    d_const_m2 = 20
    d_people_m2 = 30
    L_total_m2 = 115
    
    # Distances of noise sources from the people's location
    d_train_to_people = 30
    d_const_to_people = 50

    # --- Step 2: Calculate the source strength of each source ---
    # We work with relative intensity (I/I_0) which is 10^(L/10).
    # The 'source power' W is defined such that Relative Intensity = W / distance^2.

    # People's source power
    I_rel_people_ref = 10**(L_people_ref / 10)
    W_people = I_rel_people_ref * (d_people_ref**2)

    # Train's source power (from measurement 1)
    I_rel_total_m1 = 10**(L_total_m1 / 10)
    I_rel_people_m1 = W_people / (d_people_m1**2)
    I_rel_train_m1 = I_rel_total_m1 - I_rel_people_m1
    W_train = I_rel_train_m1 * (d_train_m1**2)

    # Construction's source power (from measurement 2)
    I_rel_total_m2 = 10**(L_total_m2 / 10)
    I_rel_people_m2 = W_people / (d_people_m2**2)
    I_rel_const_m2 = I_rel_total_m2 - I_rel_people_m2
    W_const = I_rel_const_m2 * (d_const_m2**2)

    # --- Step 3: Calculate total noise intensity at the people's location ---
    I_rel_train_at_people = W_train / (d_train_to_people**2)
    I_rel_const_at_people = W_const / (d_const_to_people**2)

    # Total noise intensity is the sum of individual noise intensities
    I_rel_noise_total = I_rel_train_at_people + I_rel_const_at_people
    
    # --- Step 4: Convert total noise intensity back to dB ---
    L_noise = 10 * math.log10(I_rel_noise_total)

    # --- Step 5: Calculate and print the SNR ---
    SNR = L_signal - L_noise
    
    print("The final Signal-to-Noise Ratio (SNR) is calculated as:")
    print(f"SNR = L_signal - L_noise")
    print(f"SNR = {L_signal} dB - {L_noise:.2f} dB = {SNR:.2f} dB")

calculate_snr()