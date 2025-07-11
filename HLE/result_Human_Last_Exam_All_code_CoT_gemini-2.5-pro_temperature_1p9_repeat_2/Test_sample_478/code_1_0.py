import math

def solve_snr():
    """
    This function calculates the Signal-to-Noise Ratio (SNR) based on the given parameters.
    """
    # Step 1: Define initial parameters from the problem description.
    # The signal is the sound level of the people talking at their location.
    L_signal = 75.0  # dB

    # Data for Measurement 1 (Train and People)
    L_total_1 = 100.0  # dB
    r_train_1 = 10.0   # meters
    r_people_1 = 20.0  # meters

    # Data for Measurement 2 (Construction and People)
    L_total_2 = 115.0  # dB
    r_const_2 = 20.0   # meters
    r_people_2 = 30.0  # meters

    # Distances of noise sources from the final location (where the people are)
    r_train_final = 30.0   # meters
    r_const_final = 50.0   # meters

    # Step 2: We assume the 75 dB level of the people acts like a source of 75 dB at 1m.
    # Let's calculate the people's sound level at the two measurement points to confirm it's negligible.
    L_source_people = 75.0 # dB at 1m

    # People's sound level at measurement point 1. The effect is negligible on 100 dB.
    # L_people_at_1 = L_source_people - 20 * math.log10(r_people_1) -> 75 - 26.02 = 48.98 dB
    # Adding 48.98 dB to a 100 dB source results in 100.00003 dB.
    # So we can approximate the train level at 10m as 100 dB.
    L_train_at_10m = 100.0

    # Similarly, the people's sound level at measurement point 2 is negligible for 115 dB.
    # L_people_at_2 = L_source_people - 20 * math.log10(r_people_2) -> 75 - 29.54 = 45.46 dB
    # Adding 45.46 dB to a 115 dB source results in 115.00001 dB.
    # So we can approximate the construction level at 20m as 115 dB.
    L_const_at_20m = 115.0

    # Step 3: Calculate the noise levels of the train and construction at the final location.
    # L2 = L1 - 20 * log10(r2 / r1)
    L_train_final = L_train_at_10m - 20 * math.log10(r_train_final / r_train_1)
    L_const_final = L_const_at_20m - 20 * math.log10(r_const_final / r_const_2)

    # Step 4: Calculate the total noise level by adding the intensities of the two noise sources.
    # I = 10^(L/10)
    I_train_final = 10**(L_train_final / 10)
    I_const_final = 10**(L_const_final / 10)
    
    # Total noise intensity is the sum of individual intensities.
    I_total_noise = I_train_final + I_const_final

    # Convert total intensity back to decibels. L = 10 * log10(I)
    L_total_noise = 10 * math.log10(I_total_noise)

    # Step 5: Calculate the Signal-to-Noise Ratio (SNR).
    SNR = L_signal - L_total_noise

    # Print the results, including the final equation as requested.
    print(f"The signal level at the location is: {L_signal:.2f} dB")
    print(f"The train noise level at the location is: {L_train_final:.2f} dB")
    print(f"The construction noise level at the location is: {L_const_final:.2f} dB")
    print(f"The total combined noise level at the location is: {L_total_noise:.2f} dB\n")
    print("The Signal-to-Noise Ratio (SNR) is calculated as Signal Level - Total Noise Level.")
    print(f"SNR = {L_signal:.2f} dB - {L_total_noise:.2f} dB = {SNR:.2f} dB")

solve_snr()