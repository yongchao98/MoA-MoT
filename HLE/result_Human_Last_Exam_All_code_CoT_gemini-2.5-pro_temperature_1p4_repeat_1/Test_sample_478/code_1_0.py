import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio at a specific location based on given sound levels and distances.
    """
    # Step 1: Define given values
    # Signal level (S) at the people's location is 75 dB.
    S_dB = 75
    # The source level of the people at a reference distance of 1m is also taken as 75 dB.
    Lp_source_1m = 75

    # Measurement 1: Train and People
    r_train_m1 = 10  # distance from train in meters
    r_people_m1 = 20 # distance from people in meters
    L_total_m1 = 100 # total sound level in dB

    # Measurement 2: Construction and People
    r_construction_m2 = 20 # distance from construction in meters
    r_people_m2 = 30   # distance from people in meters
    L_total_m2 = 115   # total sound level in dB

    # Distances of noise sources from the people's location (Point P)
    r_train_at_P = 30
    r_construction_at_P = 50

    # Step 2: Calculate the source level of the train
    # First, calculate the people's sound contribution at measurement point 1.
    Lp_at_m1 = Lp_source_1m - 20 * math.log10(r_people_m1)
    # The contribution from people (approx 49 dB) is negligible compared to 100 dB.
    # So, we can approximate the train's sound level at 10m as 100 dB.
    Lt_at_m1 = 100
    # Calculate the train's source level at 1m.
    Lt_source_1m = Lt_at_m1 + 20 * math.log10(r_train_m1)

    # Step 3: Calculate the source level of the construction site
    # First, calculate the people's sound contribution at measurement point 2.
    Lp_at_m2 = Lp_source_1m - 20 * math.log10(r_people_m2)
    # The contribution from people (approx 45 dB) is negligible compared to 115 dB.
    # So, we can approximate the construction's sound level at 20m as 115 dB.
    Lc_at_m2 = 115
    # Calculate the construction's source level at 1m.
    Lc_source_1m = Lc_at_m2 + 20 * math.log10(r_construction_m2)

    # Step 4: Calculate noise levels at the people's location (Point P)
    # Noise from the train at Point P (30m away)
    Lt_at_P = Lt_source_1m - 20 * math.log10(r_train_at_P)
    # Noise from the construction at Point P (50m away)
    Lc_at_P = Lc_source_1m - 20 * math.log10(r_construction_at_P)

    # Step 5: Calculate total noise level (N) at Point P by combining train and construction noise
    # Convert dB to intensity-like values, sum them, and convert back to dB.
    intensity_t = 10**(Lt_at_P / 10)
    intensity_c = 10**(Lc_at_P / 10)
    total_intensity = intensity_t + intensity_c
    N_dB = 10 * math.log10(total_intensity)

    # Step 6: Calculate the Signal-to-Noise Ratio (SNR)
    SNR = S_dB - N_dB
    
    # Print the final equation and result
    print(f"Signal (S) at the location = {S_dB} dB")
    print(f"Noise from Train at the location = {Lt_at_P:.2f} dB")
    print(f"Noise from Construction at the location = {Lc_at_P:.2f} dB")
    print(f"Total Noise (N) at the location = {N_dB:.2f} dB")
    print("\n--- Final Calculation ---")
    print(f"SNR = S - N")
    print(f"SNR = {S_dB} - {N_dB:.2f} = {SNR:.2f} dB")


solve_snr()
<<<D>>>