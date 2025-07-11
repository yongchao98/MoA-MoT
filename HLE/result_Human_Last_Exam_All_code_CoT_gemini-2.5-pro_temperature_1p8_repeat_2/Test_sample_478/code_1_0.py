import math

def calculate_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # 1. Define Signal
    signal_db = 75.0
    
    # Information about the noise sources
    # Train
    train_level_1 = 100.0  # dB
    train_dist_1 = 10.0    # meters
    train_dist_2 = 30.0    # meters (distance to people)

    # Construction
    const_level_1 = 115.0  # dB
    const_dist_1 = 20.0    # meters
    const_dist_2 = 50.0    # meters (distance to people)
    
    # 2. Calculate noise level from each source at the people's location
    # Formula: L2 = L1 - 20 * log10(r2/r1)
    
    # Noise from train at the people's location
    noise_train_db = train_level_1 - 20 * math.log10(train_dist_2 / train_dist_1)
    
    # Noise from construction at the people's location
    noise_const_db = const_level_1 - 20 * math.log10(const_dist_2 / const_dist_1)

    # 3. Combine the noise sources
    # Convert dB to a value proportional to intensity (power)
    # I_prop = 10^(L_db / 10)
    intensity_train = 10**(noise_train_db / 10)
    intensity_const = 10**(noise_const_db / 10)
    
    # Add intensities to get total noise intensity
    total_intensity_noise = intensity_train + intensity_const
    
    # Convert total intensity back to dB
    # L_db = 10 * log10(I)
    total_noise_db = 10 * math.log10(total_intensity_noise)

    # 4. Calculate Signal-to-Noise Ratio (SNR)
    # SNR = Signal (dB) - Noise (dB)
    snr = signal_db - total_noise_db

    # Print the steps as per the instructions
    print("Step 1: The signal level is defined as 75 dB.")
    print(f"Signal = {signal_db:.2f} dB\n")
    
    print("Step 2: Calculate the noise level from each source at the people's location.")
    print(f"Noise from train at {train_dist_2}m = {train_level_1:.2f} - 20*log10({train_dist_2}/{train_dist_1}) = {noise_train_db:.2f} dB")
    print(f"Noise from construction at {const_dist_2}m = {const_level_1:.2f} - 20*log10({const_dist_2}/{const_dist_1}) = {noise_const_db:.2f} dB\n")

    print("Step 3: Combine noise sources to find the total noise level.")
    print(f"Total Noise = 10*log10(10^({noise_train_db:.2f}/10) + 10^({noise_const_db:.2f}/10)) = {total_noise_db:.2f} dB\n")

    print("Step 4: Calculate the final Signal-to-Noise Ratio (SNR).")
    print("Final Equation: SNR = Signal_dB - Total_Noise_dB")
    print(f"SNR = {signal_db:.2f} - {total_noise_db:.2f}")
    print(f"Final SNR = {snr:.2f} dB")

# Run the calculation and print the result
calculate_snr()