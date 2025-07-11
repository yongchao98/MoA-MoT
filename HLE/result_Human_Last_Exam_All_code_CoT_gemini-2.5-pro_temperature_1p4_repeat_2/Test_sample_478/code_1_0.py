import math

def calculate_snr():
    """
    Calculates the signal-to-noise ratio at the location of a group of people.
    """
    # --- Step 1: Determine effective noise source strengths ---
    # We assume the people's sound contribution is negligible in the mixed measurements
    # due to the large difference in decibels, so we approximate.
    # L_total = 10*log10(10^(L_noise/10) + 10^(L_people/10)).
    # If L_noise is much larger than L_people, L_total is approximately L_noise.
    
    # Sound level of the train measured at 10 meters is approximately 100 dB.
    L_train_at_10m = 100.0  # dB
    # Sound level of the construction measured at 20 meters is approximately 115 dB.
    L_construction_at_20m = 115.0  # dB

    # --- Step 2: Calculate Noise Levels at the Target Location ---
    # The target location is where the people are: 30m from the train and 50m from the construction.
    # Use the distance attenuation formula: L2 = L1 - 20 * log10(r2 / r1)

    # Calculate train's noise level at 30 meters
    d_train_initial = 10.0  # meters
    d_train_final = 30.0    # meters
    L_train_at_people = L_train_at_10m - 20 * math.log10(d_train_final / d_train_initial)

    # Calculate construction's noise level at 50 meters
    d_construction_initial = 20.0  # meters
    d_construction_final = 50.0    # meters
    L_construction_at_people = L_construction_at_20m - 20 * math.log10(d_construction_final / d_construction_initial)

    # --- Step 3: Calculate Total Noise Level (L_N) ---
    # Combine the two noise sources by adding their intensities.
    # L_N = 10 * log10(10^(L_train/10) + 10^(L_construction/10))
    
    intensity_train = 10**(L_train_at_people / 10)
    intensity_construction = 10**(L_construction_at_people / 10)
    L_N = 10 * math.log10(intensity_train + intensity_construction)

    # --- Step 4: Calculate Signal-to-Noise Ratio (SNR) ---
    # The signal level (L_S) is the sound of the people talking.
    L_S = 75.0  # dB
    
    # SNR = L_S - L_N
    SNR = L_S - L_N
    
    print("Signal and Noise Calculation:")
    print(f"Signal Level (L_S): {L_S:.2f} dB")
    print(f"Noise from Train at Location: {L_train_at_people:.2f} dB")
    print(f"Noise from Construction at Location: {L_construction_at_people:.2f} dB")
    print(f"Total Noise Level (L_N): {L_N:.2f} dB")
    print("\nFinal Signal-to-Noise Ratio (SNR) Equation:")
    print(f"SNR = L_S - L_N")
    print(f"SNR = {L_S:.2f} dB - {L_N:.2f} dB")
    print(f"SNR = {SNR:.2f} dB")

calculate_snr()
<<<D>>>