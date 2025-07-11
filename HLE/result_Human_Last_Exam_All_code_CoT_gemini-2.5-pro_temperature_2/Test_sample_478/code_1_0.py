import math

def solve_snr():
    """
    This function calculates the Signal-to-Noise Ratio (SNR) based on the problem description.
    """
    # Step 1: Define known values from the problem statement.
    S = 75.0  # Signal level in dB.
    
    # Measurement 1 data
    L_total1 = 100.0
    r_train1 = 10.0
    r_people1 = 20.0
    
    # Measurement 2 data
    L_total2 = 115.0
    r_constr2 = 20.0
    r_people2 = 30.0
    
    # Distances of noise sources from the people's location (P)
    r_train_p = 30.0
    r_constr_p = 50.0

    print("Step 1: Determine the source sound levels for the noise sources.")
    # Assumption: The people's source level at 1m is 75 dB.
    L_people_1m = 75.0
    print(f"Assuming people's source sound level at 1m is {L_people_1m:.2f} dB.")

    # Step 2: Calculate the source level of the train.
    # We work with relative "power" (intensity), which is proportional to 10^(L/10).
    power_total1 = 10**(L_total1 / 10.0)
    # Calculate people's sound level and power at measurement point 1.
    L_people_at_m1 = L_people_1m - 20 * math.log10(r_people1)
    power_people_at_m1 = 10**(L_people_at_m1 / 10.0)
    # The train's power is the total minus the people's contribution.
    power_train_at_m1 = power_total1 - power_people_at_m1
    L_train_at_m1 = 10 * math.log10(power_train_at_m1)
    # Calculate train's source level at 1m.
    L_train_1m = L_train_at_m1 + 20 * math.log10(r_train1)
    print(f"Calculated train source sound level at 1m: {L_train_1m:.2f} dB.")

    # Step 3: Calculate the source level of the construction site.
    power_total2 = 10**(L_total2 / 10.0)
    # Calculate people's sound level and power at measurement point 2.
    L_people_at_m2 = L_people_1m - 20 * math.log10(r_people2)
    power_people_at_m2 = 10**(L_people_at_m2 / 10.0)
    # The construction's power is the total minus the people's contribution.
    power_constr_at_m2 = power_total2 - power_people_at_m2
    L_constr_at_m2 = 10 * math.log10(power_constr_at_m2)
    # Calculate construction's source level at 1m.
    L_constr_1m = L_constr_at_m2 + 20 * math.log10(r_constr2)
    print(f"Calculated construction source sound level at 1m: {L_constr_1m:.2f} dB.")

    print("\nStep 2: Calculate the total noise level at the people's location.")
    # Step 4: Calculate the noise level from each source at the people's location.
    L_train_at_p = L_train_1m - 20 * math.log10(r_train_p)
    L_constr_at_p = L_constr_1m - 20 * math.log10(r_constr_p)
    print(f"Noise level from train at people's location ({r_train_p}m away): {L_train_at_p:.2f} dB.")
    print(f"Noise level from construction at people's location ({r_constr_p}m away): {L_constr_at_p:.2f} dB.")

    # Step 5: Combine the two noise levels logarithmically to find the total noise N.
    power_train_at_p = 10**(L_train_at_p / 10.0)
    power_constr_at_p = 10**(L_constr_at_p / 10.0)
    total_noise_power_at_p = power_train_at_p + power_constr_at_p
    N = 10 * math.log10(total_noise_power_at_p)
    print(f"Total noise level (N) at people's location: {N:.2f} dB.")

    print("\nStep 3: Calculate the Signal-to-Noise Ratio (SNR).")
    # Step 6: Calculate the final SNR.
    SNR = S - N
    print("\nThe final equation is SNR = S - N")
    print(f"SNR = {S:.2f} dB - {N:.2f} dB = {SNR:.2f} dB")
    
if __name__ == '__main__':
    solve_snr()