import math

def calculate_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define given values and constants.
    # We interpret "people talking at 75 dB" as their source strength at 1 meter.
    L_p_source = 75.0  # dB at 1m for people
    S_db = 75.0 # Signal level in dB is given as the people's talking level.

    # Measurement 1: Total sound level of 100 dB at a location
    # 10m from the train and 20m from the people.
    L_total1 = 100.0
    d_train1 = 10.0
    d_people1 = 20.0

    # Measurement 2: Total sound level of 115 dB at a location
    # 20m from the construction and 30m from the people.
    L_total2 = 115.0
    d_construction2 = 20.0
    d_people2 = 30.0

    # Distances for the final noise calculation at the people's location
    d_train_to_people = 30.0
    d_construction_to_people = 50.0

    # It's easier to work with intensities. Intensity I is proportional to 10^(L/10).
    # Intensity from a source follows the inverse square law: I(r) = I_source / r^2
    # So, L(r) = L_source - 20 * log10(r)

    # Step 2: Use measurements to find the source strength of the train and construction.

    # From Measurement 1, find train source strength (L_t_source)
    # L_people_at_m1 = L_p_source - 20 * log10(d_people1)
    # 10^(L_total1/10) = 10^(L_train_at_m1/10) + 10^(L_people_at_m1/10)
    
    # Let's convert to relative intensity first
    I_total1_rel = 10**(L_total1 / 10)
    I_people_source_rel = 10**(L_p_source / 10)
    
    I_people_at_m1_rel = I_people_source_rel / (d_people1**2)
    I_train_at_m1_rel = I_total1_rel - I_people_at_m1_rel
    I_train_source_rel = I_train_at_m1_rel * (d_train1**2)
    L_t_source = 10 * math.log10(I_train_source_rel)

    # From Measurement 2, find construction source strength (L_c_source)
    I_total2_rel = 10**(L_total2 / 10)

    I_people_at_m2_rel = I_people_source_rel / (d_people2**2)
    I_construction_at_m2_rel = I_total2_rel - I_people_at_m2_rel
    I_construction_source_rel = I_construction_at_m2_rel * (d_construction2**2)
    L_c_source = 10 * math.log10(I_construction_source_rel)

    # Step 3: Calculate the total noise level (N_db) at the people's location.
    # Noise is the combination of the train and construction sound.
    
    # Calculate relative intensity of train noise at the people's location
    I_train_at_people_rel = I_train_source_rel / (d_train_to_people**2)

    # Calculate relative intensity of construction noise at the people's location
    I_construction_at_people_rel = I_construction_source_rel / (d_construction_to_people**2)

    # Total noise intensity is the sum of the two
    I_noise_total_rel = I_train_at_people_rel + I_construction_at_people_rel
    
    # Convert total noise intensity back to dB
    N_db = 10 * math.log10(I_noise_total_rel)

    # Step 4: Calculate the Signal-to-Noise Ratio (SNR).
    SNR_db = S_db - N_db
    
    # Print the final equation and result
    print("The final equation is: SNR = S - N")
    print(f"Signal (S) = {S_db:.2f} dB")
    print(f"Noise (N) = {N_db:.2f} dB")
    print(f"Signal-to-Noise Ratio (SNR) = {S_db:.2f} - {N_db:.2f} = {SNR_db:.2f} dB")

calculate_snr()
<<<D>>>