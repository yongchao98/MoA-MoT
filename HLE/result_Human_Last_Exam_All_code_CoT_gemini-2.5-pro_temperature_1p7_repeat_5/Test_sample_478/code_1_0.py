import math

def solve_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define initial values from the problem
    # Signal: Sound level of people talking at their location (assumed at r=1m for propagation)
    s_people_talk = 75  # dB

    # Distances of noise sources from the people's location
    r_train_to_people = 30  # meters
    r_construction_to_people = 50  # meters

    # Measurement Point 1 (M1) data
    l_total_m1 = 100  # dB
    r_train_to_m1 = 10  # meters
    r_people_to_m1 = 20  # meters

    # Measurement Point 2 (M2) data
    l_total_m2 = 115  # dB
    r_construction_to_m2 = 20  # meters
    r_people_to_m2 = 30  # meters
    
    print("This script calculates the Signal-to-Noise Ratio (SNR) at the location of the people talking.")
    print("The final calculation is SNR = Signal Level (dB) - Total Noise Level (dB).\n")

    # Step 2: Calculate the contribution of the people's speech at the two measurement points.
    # Formula: L_at_r = L_at_ref - 20 * log10(r / r_ref), assuming r_ref = 1m.
    l_people_at_m1 = s_people_talk - 20 * math.log10(r_people_to_m1)
    l_people_at_m2 = s_people_talk - 20 * math.log10(r_people_to_m2)

    # Step 3: Isolate the sound level of the train and construction at their respective measurement points.
    # We need to subtract the sound intensity of the people from the total measured intensity.
    # I = 10^(L/10)
    
    # For the train at M1
    i_total_m1 = 10**(l_total_m1 / 10)
    i_people_m1 = 10**(l_people_at_m1 / 10)
    i_train_m1 = i_total_m1 - i_people_m1
    l_train_at_m1 = 10 * math.log10(i_train_m1)

    # For the construction at M2
    i_total_m2 = 10**(l_total_m2 / 10)
    i_people_m2 = 10**(l_people_at_m2 / 10)
    i_construction_m2 = i_total_m2 - i_people_m2
    l_construction_at_m2 = 10 * math.log10(i_construction_m2)
    
    # Step 4: Calculate the noise level from each source at the people's location.
    
    # Noise from the train
    l_train_at_people = l_train_at_m1 - 20 * math.log10(r_train_to_people / r_train_to_m1)
    
    # Noise from the construction
    l_construction_at_people = l_construction_at_m2 - 20 * math.log10(r_construction_to_people / r_construction_to_m2)

    # Step 5: Combine the two noise levels at the people's location to find the total noise level (N).
    i_train_at_people = 10**(l_train_at_people / 10)
    i_construction_at_people = 10**(l_construction_at_people / 10)
    i_total_noise_at_people = i_train_at_people + i_construction_at_people
    n_total_noise = 10 * math.log10(i_total_noise_at_people)

    # Step 6: Calculate the Signal-to-Noise Ratio (SNR).
    snr = s_people_talk - n_total_noise
    
    print("1. The Signal (S) is the sound of the people talking, given as {} dB.".format(s_people_talk))
    
    print("\n2. The Total Noise (N) is the sum of noise from the train and construction.")
    print("   - Noise from train at people's location = {:.2f} dB".format(l_train_at_people))
    print("   - Noise from construction at people's location = {:.2f} dB".format(l_construction_at_people))
    print("   - Combining these two noise sources gives a Total Noise level of {:.2f} dB.".format(n_total_noise))
    
    print("\n3. The Signal-to-Noise Ratio (SNR) is calculated as:")
    print("   SNR = Signal - Noise")
    print("   SNR = {} - {:.2f} = {:.2f} dB".format(s_people_talk, n_total_noise, snr))

solve_snr()