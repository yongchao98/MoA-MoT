import math

def calculate_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """

    # --- Part 1: Define constants and assumptions from the problem ---
    # Assumption: "people talking at 75 dB" means their source level at 1m is 75 dB.
    # This value is also used as the signal level in the final SNR calculation.
    L_people_source_1m = 75.0  # dB

    # Measurement 1: For finding the train's source level
    L_mixed_1 = 100.0  # dB
    r_train_1 = 10.0   # meters
    r_people_1 = 20.0  # meters

    # Measurement 2: For finding the construction's source level
    L_mixed_2 = 115.0  # dB
    r_const_2 = 20.0   # meters
    r_people_2 = 30.0  # meters
    
    # Distances for final noise calculation at the people's location
    r_train_to_people = 30.0  # meters
    r_const_to_people = 50.0  # meters
    
    # --- Part 2: Calculate the source level of the train (L_T1) ---
    # Sound level of people at measurement point 1
    L_people_at_M1 = L_people_source_1m - 20 * math.log10(r_people_1)
    # Convert levels to intensities (relative to reference intensity)
    I_total_at_M1 = 10**(L_mixed_1 / 10)
    I_people_at_M1 = 10**(L_people_at_M1 / 10)
    # Find train's intensity and convert back to level
    I_train_at_M1 = I_total_at_M1 - I_people_at_M1
    L_train_at_M1 = 10 * math.log10(I_train_at_M1)
    # Calculate train's source level at 1m
    L_train_source_1m = L_train_at_M1 + 20 * math.log10(r_train_1)
    
    # --- Part 3: Calculate the source level of the construction site (L_C1) ---
    # Sound level of people at measurement point 2
    L_people_at_M2 = L_people_source_1m - 20 * math.log10(r_people_2)
    # Convert levels to intensities
    I_total_at_M2 = 10**(L_mixed_2 / 10)
    I_people_at_M2 = 10**(L_people_at_M2 / 10)
    # Find construction's intensity and convert back to level
    I_const_at_M2 = I_total_at_M2 - I_people_at_M2
    L_const_at_M2 = 10 * math.log10(I_const_at_M2)
    # Calculate construction's source level at 1m
    L_const_source_1m = L_const_at_M2 + 20 * math.log10(r_const_2)
    
    # --- Part 4: Calculate the total noise level at the people's location ---
    # Noise level from train at the people's location
    L_train_at_people = L_train_source_1m - 20 * math.log10(r_train_to_people)
    # Noise level from construction at the people's location
    L_const_at_people = L_const_source_1m - 20 * math.log10(r_const_to_people)
    # Convert individual noise levels to intensities
    I_train_noise = 10**(L_train_at_people / 10)
    I_const_noise = 10**(L_const_at_people / 10)
    # Add intensities and convert back to a total noise level in dB
    L_total_noise = 10 * math.log10(I_train_noise + I_const_noise)

    # --- Part 5: Calculate the Signal-to-Noise Ratio (SNR) ---
    L_signal = L_people_source_1m # The signal is given as 75 dB
    snr = L_signal - L_total_noise
    
    # --- Part 6: Print the results ---
    print("Calculating Signal-to-Noise Ratio (SNR)\n")
    print(f"Signal Level (S) = {L_signal:.2f} dB")
    print("Noise is the sum of sound from the train and construction.")
    print(f"Noise from train at people's location = {L_train_at_people:.2f} dB")
    print(f"Noise from construction at people's location = {L_const_at_people:.2f} dB")
    print("\nTo find Total Noise Level (N), we add their intensities:")
    print(f"N = 10 * log10( 10^({L_train_at_people:.2f}/10) + 10^({L_const_at_people:.2f}/10) )")
    print(f"N = {L_total_noise:.2f} dB")
    print("\nFinal SNR Calculation:")
    print(f"SNR = S - N")
    print(f"SNR = {L_signal:.2f} dB - {L_total_noise:.2f} dB = {snr:.2f} dB")


if __name__ == '__main__':
    calculate_snr()
<<<D>>>