import math

def calculate_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define given values from the problem statement
    # Signal (S) is the sound level of the people at their own location.
    S_people_signal = 75.0  # dB

    # Measurement point 1 data
    L_total1 = 100.0  # dB
    r_train1 = 10.0   # meters
    r_people1 = 20.0  # meters

    # Measurement point 2 data
    L_total2 = 115.0  # dB
    r_const2 = 20.0   # meters
    r_people2 = 30.0  # meters

    # Distances from noise sources to the people's location
    r_train_to_people = 30.0 # meters
    r_const_to_people = 50.0 # meters
    
    # Assumption: The people's source sound level at r=1m is 75 dB.
    L_people_source = 75.0 # dB at 1m

    # Step 2: Calculate the Train's source level
    # First, find the sound level contributed by the people at measurement point 1
    # Sound Propagation Formula: L(r) = L_source - 20 * log10(r)
    L_people_at_p1 = L_people_source - 20 * math.log10(r_people1)
    
    # Convert decibels to intensities to perform subtraction
    # Intensity formula: I = 10^(L/10)
    I_total1 = 10**(L_total1 / 10)
    I_people_at_p1 = 10**(L_people_at_p1 / 10)
    
    # The train's intensity is the total minus the people's
    I_train_at_p1 = I_total1 - I_people_at_p1
    
    # Convert train's intensity back to dB
    # Decibel formula: L = 10 * log10(I)
    L_train_at_p1 = 10 * math.log10(I_train_at_p1)
    
    # Now, calculate the train's source level at 1 meter
    L_train_source = L_train_at_p1 + 20 * math.log10(r_train1)

    # Step 3: Calculate the Construction's source level (similar process)
    L_people_at_p2 = L_people_source - 20 * math.log10(r_people2)
    I_total2 = 10**(L_total2 / 10)
    I_people_at_p2 = 10**(L_people_at_p2 / 10)
    I_const_at_p2 = I_total2 - I_people_at_p2
    L_const_at_p2 = 10 * math.log10(I_const_at_p2)
    L_const_source = L_const_at_p2 + 20 * math.log10(r_const2)

    # Step 4: Calculate total noise level at the people's location
    # Noise from train at the people's location (distance = 30m)
    N_train = L_train_source - 20 * math.log10(r_train_to_people)
    
    # Noise from construction at the people's location (distance = 50m)
    N_const = L_const_source - 20 * math.log10(r_const_to_people)
    
    # Combine the two noises. We must add their intensities, not their dB values.
    I_noise_train = 10**(N_train / 10)
    I_noise_const = 10**(N_const / 10)
    I_total_noise = I_noise_train + I_noise_const
    N_total = 10 * math.log10(I_total_noise)
    
    # Step 5: Calculate the final Signal-to-Noise Ratio (SNR)
    SNR = S_people_signal - N_total
    
    # Print the results step-by-step
    print(f"1. People's Source Level (assumed at 1m): {L_people_source:.2f} dB")
    print(f"2. Train's Source Level (at 1m): {L_train_source:.2f} dB")
    print(f"3. Construction's Source Level (at 1m): {L_const_source:.2f} dB")
    print("\n--- At the people's location ---")
    print(f"Noise from Train (at {r_train_to_people}m): {N_train:.2f} dB")
    print(f"Noise from Construction (at {r_const_to_people}m): {N_const:.2f} dB")
    print(f"Total Combined Noise Level: {N_total:.2f} dB")
    print("\n--- Final Signal-to-Noise Ratio (SNR) ---")
    print(f"Signal = {S_people_signal:.2f} dB")
    print(f"Noise = {N_total:.2f} dB")
    print(f"SNR = Signal - Noise")
    print(f"SNR = {S_people_signal:.2f} - {N_total:.2f} = {SNR:.2f} dB")

calculate_snr()