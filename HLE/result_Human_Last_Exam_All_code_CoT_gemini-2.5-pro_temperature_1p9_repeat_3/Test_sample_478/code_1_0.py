import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define given variables from the problem.
    # The signal level (L_S) from the people at the point of interest (POI).
    L_S = 75.0  # dB

    # Distances of noise sources to the POI.
    d_train_poi = 30.0  # meters
    d_const_poi = 50.0  # meters

    # Measurement 1 parameters.
    L_m1 = 100.0  # dB
    d_train_m1 = 10.0  # meters
    d_people_m1 = 20.0  # meters

    # Measurement 2 parameters.
    L_m2 = 115.0  # dB
    d_const_m2 = 20.0  # meters
    d_people_m2 = 30.0  # meters

    print("Step 1: Determine the source strength of each sound source.")
    # We assume "people talking at 75 dB" means their sound level at a 1-meter distance is 75 dB.
    # From this, we can calculate their source strength, S_people.
    # I_prime = 10^(L/10), S = I_prime * r^2
    L_people_ref = 75.0
    r_ref = 1.0
    S_people = (10**(L_people_ref / 10)) * (r_ref**2)
    print(f"  - Source strength of people (S_people) based on 75 dB at 1m: {S_people:.2e}")

    # Now solve for the other two source strengths using the two measurements.
    # Measurement 1: S_train/d_train_m1^2 + S_people/d_people_m1^2 = 10^(L_m1/10)
    I_prime_total_m1 = 10**(L_m1 / 10)
    S_train = (I_prime_total_m1 - S_people / (d_people_m1**2)) * (d_train_m1**2)
    print(f"  - Source strength of train (S_train): {S_train:.2e}")
    
    # Measurement 2: S_construction/d_const_m2^2 + S_people/d_people_m2^2 = 10^(L_m2/10)
    I_prime_total_m2 = 10**(L_m2 / 10)
    S_construction = (I_prime_total_m2 - S_people / (d_people_m2**2)) * (d_const_m2**2)
    print(f"  - Source strength of construction (S_construction): {S_construction:.2e}\n")

    print("Step 2: Calculate the total noise intensity at the point of interest (POI).")
    # Noise at POI is the sum of intensities from the train (30m away) and construction (50m away).
    I_prime_train_poi = S_train / (d_train_poi**2)
    I_prime_const_poi = S_construction / (d_const_poi**2)
    I_prime_noise_total = I_prime_train_poi + I_prime_const_poi
    print(f"  - Combined relative noise intensity at POI: {I_prime_noise_total:.2e}\n")
    
    print("Step 3: Convert total noise intensity to a sound level (L_N) in decibels.")
    L_N = 10 * math.log10(I_prime_noise_total)
    print(f"  - Total noise level (L_N): {L_N:.2f} dB\n")

    print("Step 4: Calculate the final Signal-to-Noise Ratio (SNR).")
    # SNR (dB) = Signal Level (dB) - Noise Level (dB)
    SNR = L_S - L_N
    
    # Print the final equation with all numbers
    print(f"The final calculation is:")
    print(f"SNR = L_S - L_N")
    print(f"SNR = {L_S} dB - {L_N:.2f} dB = {SNR:.2f} dB")


solve_snr()