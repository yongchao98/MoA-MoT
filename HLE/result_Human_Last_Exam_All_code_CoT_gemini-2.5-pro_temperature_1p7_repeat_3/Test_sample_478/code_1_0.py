import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio at the location of the people.
    """
    # Step 1: Define the given values from the problem.
    # The signal is the sound level of the people at their location.
    L_signal = 75  # dB

    # Measurement 1: Train and People
    L_total_1 = 100  # dB
    r_train_1 = 10   # meters from train
    r_people_1 = 20  # meters from people

    # Measurement 2: Construction and People
    L_total_2 = 115  # dB
    r_construction_2 = 20 # meters from construction
    r_people_2 = 30       # meters from people
    
    # We assume the 75 dB of the people is the source level at a reference distance of 1m.
    # This is a standard assumption when a source level at the source isn't specified.
    L_people_ref = 75 # dB
    r_people_ref = 1  # m

    # Location of interest (where the people are)
    r_train_final = 30      # meters from train
    r_construction_final = 50 # meters from construction

    # Step 2: Determine the sound level of each noise source from the measurements.
    # For measurement 1, let's find the sound level of just the train at 10m.
    # First, calculate the people's contribution at this point.
    L_people_at_m1 = L_people_ref - 20 * math.log10(r_people_1 / r_people_ref)
    
    # We can see the people's contribution is negligible compared to 100 dB,
    # so we can approximate the train's sound level at 10m as 100 dB.
    L_train_at_10m = L_total_1

    # For measurement 2, let's find the sound level of the construction site at 20m.
    # The people's contribution is also negligible here compared to 115 dB.
    L_construction_at_20m = L_total_2
    
    # Step 3: Calculate the sound level of each noise source at the people's location.
    # Sound level of the train at 30 meters.
    L_train_final = L_train_at_10m - 20 * math.log10(r_train_final / r_train_1)

    # Sound level of the construction site at 50 meters.
    L_construction_final = L_construction_at_20m - 20 * math.log10(r_construction_final / r_construction_2)

    # Step 4: Combine the noise levels to get the total noise level at the location.
    # We must convert dB to intensity, add them, and then convert back to dB.
    I_train_final = 10**(L_train_final / 10)
    I_construction_final = 10**(L_construction_final / 10)
    
    I_total_noise = I_train_final + I_construction_final
    L_total_noise = 10 * math.log10(I_total_noise)

    # Step 5: Calculate the Signal-to-Noise Ratio (SNR).
    SNR = L_signal - L_total_noise

    # Print the results step-by-step
    print(f"Signal Level (L_Signal): {L_signal} dB")
    print("-" * 30)
    print(f"Noise from train at the location ({r_train_final}m):")
    print(f"L_train = {L_train_at_10m} dB - 20*log10({r_train_final}/{r_train_1}) = {L_train_final:.2f} dB")
    print(f"Noise from construction at the location ({r_construction_final}m):")
    print(f"L_construction = {L_construction_at_20m} dB - 20*log10({r_construction_final}/{r_construction_2}) = {L_construction_final:.2f} dB")
    print("-" * 30)
    print("Combining noise sources:")
    print(f"Total Noise Level (L_Noise) = 10*log10(10^({L_train_final:.2f}/10) + 10^({L_construction_final:.2f}/10)) = {L_total_noise:.2f} dB")
    print("-" * 30)
    print("Final Signal-to-Noise Ratio (SNR) Calculation:")
    print(f"SNR = L_Signal - L_Noise")
    print(f"SNR = {L_signal} dB - {L_total_noise:.2f} dB = {SNR:.2f} dB")

solve_snr()