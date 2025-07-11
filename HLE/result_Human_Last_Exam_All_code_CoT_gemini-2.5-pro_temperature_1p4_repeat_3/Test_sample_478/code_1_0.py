import math

def calculate_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define known values from the problem statement
    S = 75  # Signal level in dB at the people's location

    # Distances for the final noise calculation
    r_train_to_people = 30
    r_construction_to_people = 50

    # Data Point A for train calculation
    L_A = 100
    r_train_A = 10
    r_people_A = 20

    # Data Point B for construction calculation
    L_B = 115
    r_construction_B = 20
    r_people_B = 30

    # Step 2: Assumption: The source sound level of the people at 1m is 75 dB
    Lp_1m = 75
    print(f"Step 1: The signal level (S) is given as {S} dB.")
    print("Assuming the source sound level of the people is 75 dB at 1 meter.\n")

    # Step 3: Calculate the source sound level of the train
    print("Step 2: Calculate the source level of the train.")
    # The contribution from the people at measurement point A is negligible,
    # so we can approximate the train's sound level at 10m as 100 dB.
    # L_train(10m) ≈ 100 dB
    # A more precise calculation:
    Lp_A = Lp_1m - 20 * math.log10(r_people_A)
    Ip_A = 10**(Lp_A / 10)
    Itotal_A = 10**(L_A / 10)
    It_A = Itotal_A - Ip_A
    Lt_A = 10 * math.log10(It_A)
    # Calculate the train's source level at 1 meter
    Lt_1m = Lt_A + 20 * math.log10(r_train_A)
    print(f"The calculated source sound level of the train at 1m is {Lt_1m:.2f} dB.")

    # Step 4: Calculate the source sound level of the construction
    print("\nStep 3: Calculate the source level of the construction.")
    # The contribution from the people at measurement point B is negligible,
    # so we can approximate the construction's sound level at 20m as 115 dB.
    # L_construction(20m) ≈ 115 dB
    # A more precise calculation:
    Lp_B = Lp_1m - 20 * math.log10(r_people_B)
    Ip_B = 10**(Lp_B / 10)
    Itotal_B = 10**(L_B / 10)
    Ic_B = Itotal_B - Ip_B
    Lc_B = 10 * math.log10(Ic_B)
    # Calculate the construction's source level at 1 meter
    Lc_1m = Lc_B + 20 * math.log10(r_construction_B)
    print(f"The calculated source sound level of the construction at 1m is {Lc_1m:.2f} dB.")

    # Step 5: Calculate noise levels at the people's location
    print("\nStep 4: Calculate the noise level from each source at the people's location.")
    L_train_at_people = Lt_1m - 20 * math.log10(r_train_to_people)
    L_constr_at_people = Lc_1m - 20 * math.log10(r_construction_to_people)
    print(f"Noise from train (at {r_train_to_people}m) = {L_train_at_people:.2f} dB")
    print(f"Noise from construction (at {r_construction_to_people}m) = {L_constr_at_people:.2f} dB")

    # Step 6: Combine the noise levels
    print("\nStep 5: Combine the noise levels to find the total noise (N).")
    I_train_noise = 10**(L_train_at_people / 10)
    I_constr_noise = 10**(L_constr_at_people / 10)
    N = 10 * math.log10(I_train_noise + I_constr_noise)
    print(f"Total noise level (N) = 10 * log10(10^({L_train_at_people:.2f}/10) + 10^({L_constr_at_people:.2f}/10)) = {N:.2f} dB")

    # Step 7: Calculate the final SNR
    print("\nStep 6: Calculate the Signal-to-Noise Ratio (SNR).")
    SNR = S - N
    print(f"Final Equation: SNR = S - N")
    print(f"SNR = {S} - {N:.2f} = {SNR:.2f} dB")
    return SNR

# Run the calculation and store the result
final_snr = calculate_snr()
print(f"\nThe final Signal-to-Noise ratio is approximately {final_snr:.2f} dB.")
print("<<<D>>>")
