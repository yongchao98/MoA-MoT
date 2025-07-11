import math

def calculate_snr():
    """
    This script calculates the signal-to-noise ratio based on the provided sound level measurements.
    """
    # Step 1: Define initial knowns and assumptions
    S_dB = 75.0  # Signal level in dB.
    # We assume the 75 dB value for the people is their source strength at 1 meter.
    L_people_1m = 75.0
    
    print("Step 1: Define Signal and Source Strengths")
    print(f"Signal (S) at the location = {S_dB} dB")
    print(f"Assumed source strength of people at 1m = {L_people_1m} dB\n")

    # Step 2: Use Measurement 1 to find Train's source level
    L_total_1 = 100.0
    d_train_1 = 10.0
    d_people_1 = 20.0

    # Calculate people's sound level at measurement point 1
    L_people_at_p1 = L_people_1m - 20 * math.log10(d_people_1)
    
    # Convert levels to normalized intensities
    I_total_1_norm = 10**(L_total_1 / 10)
    I_people_at_p1_norm = 10**(L_people_at_p1 / 10)

    # Subtract to get train's intensity at point 1
    I_train_at_p1_norm = I_total_1_norm - I_people_at_p1_norm
    L_train_at_p1 = 10 * math.log10(I_train_at_p1_norm)

    # Calculate train's source level at 1m
    L_train_1m = L_train_at_p1 + 20 * math.log10(d_train_1)
    
    print("Step 2: Calculate Train's Source Strength")
    print(f"From measurement 1 (100 dB at 10m from train, 20m from people):")
    print(f"Calculated Train source strength at 1m = {L_train_1m:.2f} dB\n")

    # Step 3: Use Measurement 2 to find Construction's source level
    L_total_2 = 115.0
    d_const_2 = 20.0
    d_people_2 = 30.0

    # Calculate people's sound level at measurement point 2
    L_people_at_p2 = L_people_1m - 20 * math.log10(d_people_2)

    # Convert levels to normalized intensities
    I_total_2_norm = 10**(L_total_2 / 10)
    I_people_at_p2_norm = 10**(L_people_at_p2 / 10)

    # Subtract to get construction's intensity at point 2
    I_const_at_p2_norm = I_total_2_norm - I_people_at_p2_norm
    L_const_at_p2 = 10 * math.log10(I_const_at_p2_norm)

    # Calculate construction's source level at 1m
    L_const_1m = L_const_at_p2 + 20 * math.log10(d_const_2)

    print("Step 3: Calculate Construction's Source Strength")
    print(f"From measurement 2 (115 dB at 20m from construction, 30m from people):")
    print(f"Calculated Construction source strength at 1m = {L_const_1m:.2f} dB\n")
    
    # Step 4: Calculate noise levels at the people's location
    d_train_final = 30.0
    d_const_final = 50.0

    # Noise from train at people's location
    L_train_noise = L_train_1m - 20 * math.log10(d_train_final)
    
    # Noise from construction at people's location
    L_const_noise = L_const_1m - 20 * math.log10(d_const_final)

    print("Step 4: Calculate Individual Noise Levels at People's Location")
    print(f"Noise from train (30m away) = {L_train_noise:.2f} dB")
    print(f"Noise from construction (50m away) = {L_const_noise:.2f} dB\n")

    # Step 5: Combine noise sources to get total noise N_dB
    I_train_final_norm = 10**(L_train_noise / 10)
    I_const_final_norm = 10**(L_const_noise / 10)

    I_noise_total_norm = I_train_final_norm + I_const_final_norm
    N_dB = 10 * math.log10(I_noise_total_norm)
    
    print("Step 5: Calculate Total Noise (N)")
    print(f"Total Noise Level (N) = {N_dB:.2f} dB\n")

    # Step 6: Calculate SNR
    SNR = S_dB - N_dB
    
    print("Step 6: Final Signal-to-Noise Ratio (SNR) Calculation")
    print("SNR = Signal - Noise")
    print(f"The final equation is: {SNR:.2f} dB = {S_dB:.2f} dB - {N_dB:.2f} dB")
    
    # Return the final value for verification against choices
    return SNR

if __name__ == "__main__":
    calculate_snr()