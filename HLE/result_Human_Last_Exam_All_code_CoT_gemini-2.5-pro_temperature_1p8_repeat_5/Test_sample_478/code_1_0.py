import math

def calculate_snr():
    """
    Calculates the signal-to-noise ratio based on the problem description.
    """
    # --- Step 1: Define given values ---
    S_dB = 75.0  # Signal level at the people's location in dB
    # Assumption: The people's source strength is 75 dB at 1m.
    Lp_source_1m = 75.0
    
    # Distances for the final calculation
    d_train_to_people = 30.0
    d_construction_to_people = 50.0
    
    # Measurement point 1 (M1)
    L_total_m1 = 100.0
    d_train_to_m1 = 10.0
    d_people_to_m1 = 20.0
    
    # Measurement point 2 (M2)
    L_total_m2 = 115.0
    d_construction_to_m2 = 20.0
    d_people_to_m2 = 30.0

    # --- Step 2: Analyze Measurement 1 to find Train noise level at the people's location ---
    
    # a. Calculate the people's sound level at M1
    Lp_at_m1 = Lp_source_1m - 20 * math.log10(d_people_to_m1 / 1.0)
    
    # b. Isolate the train's sound level at M1 by subtracting the people's sound intensity
    # I = 10^(L/10), where I is intensity and L is sound level in dB
    I_total_m1 = 10**(L_total_m1 / 10)
    Ip_at_m1 = 10**(Lp_at_m1 / 10)
    It_at_m1 = I_total_m1 - Ip_at_m1
    Lt_at_m1 = 10 * math.log10(It_at_m1)
    
    # c. Calculate the train's noise level at the people's location
    # L2 = L1 - 20 * log10(d2/d1)
    N_train = Lt_at_m1 - 20 * math.log10(d_train_to_people / d_train_to_m1)
    
    # --- Step 3: Analyze Measurement 2 to find Construction noise level at the people's location ---
    
    # a. Calculate the people's sound level at M2
    Lp_at_m2 = Lp_source_1m - 20 * math.log10(d_people_to_m2 / 1.0)
    
    # b. Isolate the construction's sound level at M2
    I_total_m2 = 10**(L_total_m2 / 10)
    Ip_at_m2 = 10**(Lp_at_m2 / 10)
    Ic_at_m2 = I_total_m2 - Ip_at_m2
    Lc_at_m2 = 10 * math.log10(Ic_at_m2)
    
    # c. Calculate the construction's noise level at the people's location
    N_construction = Lc_at_m2 - 20 * math.log10(d_construction_to_people / d_construction_to_m2)

    # --- Step 4: Combine the two noise levels at the people's location ---
    # L_total = 10 * log10(10^(L1/10) + 10^(L2/10))
    I_noise_train = 10**(N_train / 10)
    I_noise_construction = 10**(N_construction / 10)
    I_total_noise = I_noise_train + I_noise_construction
    N_total_dB = 10 * math.log10(I_total_noise)

    # --- Step 5: Calculate the Signal-to-Noise Ratio (SNR) ---
    snr = S_dB - N_total_dB

    # --- Print the results ---
    print(f"Signal Level (S): {S_dB:.2f} dB")
    print(f"Noise from Train at people's location: {N_train:.2f} dB")
    print(f"Noise from Construction at people's location: {N_construction:.2f} dB")
    print(f"Total Noise Level (N): {N_total_dB:.2f} dB")
    print("\nFinal Calculation:")
    print(f"SNR = Signal - Total Noise")
    print(f"SNR = {S_dB:.2f} dB - {N_total_dB:.2f} dB = {snr:.2f} dB")

calculate_snr()