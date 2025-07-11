import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio at the location of the people.
    """
    # Step 1: Define Signal and distances
    signal_db = 75.0  # dB, Signal level at the people's location

    # Distances for measurement 1
    r_train_m1 = 10.0  # meters
    # Distances for measurement 2
    r_const_m1 = 20.0  # meters
    
    # Distances from sources to the people's location
    r_train_to_people = 30.0  # meters
    r_const_to_people = 50.0  # meters
    
    # Combined sound levels from measurements
    l_combined_train_people = 100.0  # dB
    l_combined_const_people = 115.0  # dB
    
    print("This script calculates the Signal-to-Noise Ratio (SNR) at a specific location.")
    print(f"The signal level from the people is given as {signal_db} dB.\n")

    # Step 2: Determine sound level of each noise source at a known distance.
    # The contribution of the people's speech to the combined measurements is negligible.
    # So, we can approximate the source levels from the measurements.
    l_train_at_10m = l_combined_train_people
    l_const_at_20m = l_combined_const_people
    
    print("--- Calculating Noise from the Train ---")
    print(f"The noise level from the train at {r_train_m1}m is approximated as {l_train_at_10m} dB.")
    
    # Step 3: Calculate the train's noise level at the people's location (30m away)
    # Using the formula: L2 = L1 - 20 * log10(r2 / r1)
    l_train_at_people = l_train_at_10m - 20 * math.log10(r_train_to_people / r_train_m1)
    print(f"The noise level from the train at the people's location ({r_train_to_people}m) is:")
    print(f"L_train = {l_train_at_10m:.2f} - 20 * log10({r_train_to_people} / {r_train_m1}) = {l_train_at_people:.2f} dB\n")

    # Step 4: Calculate the construction's noise level at the people's location (50m away)
    print("--- Calculating Noise from the Construction Site ---")
    print(f"The noise level from construction at {r_const_m1}m is approximated as {l_const_at_20m} dB.")
    l_const_at_people = l_const_at_20m - 20 * math.log10(r_const_to_people / r_const_m1)
    print(f"The noise level from construction at the people's location ({r_const_to_people}m) is:")
    print(f"L_const = {l_const_at_20m:.2f} - 20 * log10({r_const_to_people} / {r_const_m1}) = {l_const_at_people:.2f} dB\n")

    # Step 5: Calculate the total noise level by combining the train and construction noise
    # Using the formula: L_total = 10 * log10(10^(L1/10) + 10^(L2/10))
    term1 = 10**(l_train_at_people / 10)
    term2 = 10**(l_const_at_people / 10)
    total_noise_db = 10 * math.log10(term1 + term2)
    print("--- Calculating Total Noise Level ---")
    print("The total noise is the combination of the train noise and construction noise.")
    print(f"Total Noise = 10 * log10(10^({l_train_at_people:.2f}/10) + 10^({l_const_at_people:.2f}/10)) = {total_noise_db:.2f} dB\n")
    
    # Step 6: Calculate the Signal-to-Noise Ratio (SNR)
    snr = signal_db - total_noise_db
    print("--- Calculating Final Signal-to-Noise Ratio (SNR) ---")
    print(f"SNR = Signal Level - Total Noise Level")
    print(f"SNR = {signal_db:.2f} dB - {total_noise_db:.2f} dB = {snr:.2f} dB")

solve_snr()
<<<D>>>