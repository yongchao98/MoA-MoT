import math

def calculate_snr():
    """
    Calculates the signal-to-noise ratio based on the problem description.
    """
    # L(r) = L_source - 20 * log10(r)
    # L_total = 10 * log10(10**(L1/10) + 10**(L2/10))

    # --- Step 1: Define knowns and assumptions ---
    # Signal level at the location of the people
    L_signal = 75  # dB
    # Assumption: The source level of the people at 1 meter is also 75 dB.
    L_people_source = 75  # dB at 1m

    # Data Point 1: Train and People
    L_total1 = 100  # dB
    r_train1 = 10   # m
    r_people1 = 20  # m

    # Data Point 2: Construction and People
    L_total2 = 115  # dB
    r_construction2 = 20 # m
    r_people2 = 30 # m
    
    # --- Step 2: Calculate Train's source level ---
    # People's level at measurement point 1
    L_people_at_p1 = L_people_source - 20 * math.log10(r_people1)
    
    # Convert levels to intensities to find the train's intensity at point 1
    I_total1 = 10**(L_total1 / 10)
    I_people_at_p1 = 10**(L_people_at_p1 / 10)
    I_train_at_p1 = I_total1 - I_people_at_p1
    
    # Convert train's intensity back to level and calculate its source level
    L_train_at_p1 = 10 * math.log10(I_train_at_p1)
    L_train_source = L_train_at_p1 + 20 * math.log10(r_train1)

    # --- Step 3: Calculate Construction's source level ---
    # People's level at measurement point 2
    L_people_at_p2 = L_people_source - 20 * math.log10(r_people2)

    # Convert levels to intensities to find the construction's intensity at point 2
    I_total2 = 10**(L_total2 / 10)
    I_people_at_p2 = 10**(L_people_at_p2 / 10)
    I_construction_at_p2 = I_total2 - I_people_at_p2

    # Convert construction's intensity back to level and calculate its source level
    L_construction_at_p2 = 10 * math.log10(I_construction_at_p2)
    L_construction_source = L_construction_at_p2 + 20 * math.log10(r_construction2)

    # --- Step 4: Calculate total noise level at the people's location ---
    # Distances of noise sources from the people's location
    r_train_final = 30 # m
    r_construction_final = 50 # m
    
    # Calculate level of each noise source at the final location
    L_train_final = L_train_source - 20 * math.log10(r_train_final)
    L_construction_final = L_construction_source - 20 * math.log10(r_construction_final)

    # Combine the noise levels
    I_train_final = 10**(L_train_final / 10)
    I_construction_final = 10**(L_construction_final / 10)
    I_noise_total = I_train_final + I_construction_final
    L_noise = 10 * math.log10(I_noise_total)

    # --- Step 5: Calculate and print the Signal-to-Noise Ratio (SNR) ---
    snr = L_signal - L_noise

    print("The final equation is: SNR = Signal Level - Noise Level")
    print(f"Signal Level = {L_signal:.2f} dB")
    print(f"Noise Level = {L_noise:.2f} dB")
    print(f"Final SNR = {L_signal:.2f} dB - {L_noise:.2f} dB = {snr:.2f} dB")

if __name__ == '__main__':
    calculate_snr()
<<<D>>>