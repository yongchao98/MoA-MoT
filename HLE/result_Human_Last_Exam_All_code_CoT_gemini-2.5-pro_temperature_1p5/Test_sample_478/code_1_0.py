import math

def solve_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define constants and given values from the problem.
    L_signal = 75  # dB, signal level at the people's location.

    # Distances of noise sources to the people's location.
    r_train_to_people = 30  # meters
    r_constr_to_people = 50  # meters

    # Data from the first measurement point (P1).
    L_total_p1 = 100  # dB
    r_train_to_p1 = 10  # meters

    # Data from the second measurement point (P2).
    L_total_p2 = 115  # dB
    r_constr_to_p2 = 20  # meters

    print("Step 1: Calculate the source level of the train at 1 meter.")
    # We assume the noise measured at P1 is dominated by the train.
    # The formula to find source level (L_src) at 1m from a measurement (L_at_r) at distance r is:
    # L_src = L_at_r + 20 * log10(r)
    L_train_src = L_total_p1 + 20 * math.log10(r_train_to_p1)
    print(f"The source level of the train is calculated to be {L_train_src:.2f} dB at 1m.\n")

    print("Step 2: Calculate the source level of the construction site at 1 meter.")
    # Similarly, we assume the noise measured at P2 is dominated by the construction.
    L_constr_src = L_total_p2 + 20 * math.log10(r_constr_to_p2)
    print(f"The source level of the construction site is calculated to be {L_constr_src:.2f} dB at 1m.\n")

    print("Step 3: Calculate the noise level from each source at the people's location.")
    # The formula to find the level at a distance r from a source is:
    # L_at_r = L_src - 20 * log10(r)
    L_train_at_people = L_train_src - 20 * math.log10(r_train_to_people)
    L_constr_at_people = L_constr_src - 20 * math.log10(r_constr_to_people)
    print(f"The noise from the train at the people's location (30m away) is {L_train_at_people:.2f} dB.")
    print(f"The noise from the construction site at the people's location (50m away) is {L_constr_at_people:.2f} dB.\n")

    print("Step 4: Combine the noise levels to find the total noise.")
    # Convert dB levels to relative intensities, add them, then convert back to dB.
    # Relative intensity i = 10^(L/10)
    i_train = 10**(L_train_at_people / 10)
    i_constr = 10**(L_constr_at_people / 10)
    i_total_noise = i_train + i_constr
    L_noise_total = 10 * math.log10(i_total_noise)
    print(f"The total combined noise level at the people's location is {L_noise_total:.2f} dB.\n")

    print("Step 5: Calculate the final Signal-to-Noise Ratio (SNR).")
    # SNR (dB) = Signal Level (dB) - Noise Level (dB)
    SNR = L_signal - L_noise_total
    
    # Print the final equation with all the numbers
    print("The final equation is:")
    print(f"SNR = Signal Level - Total Noise Level")
    print(f"SNR = {L_signal:.2f} dB - {L_noise_total:.2f} dB = {SNR:.2f} dB")


solve_snr()