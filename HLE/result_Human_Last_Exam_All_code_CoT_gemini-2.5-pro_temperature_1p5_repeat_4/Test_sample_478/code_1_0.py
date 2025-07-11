import math

def calculate_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define Signal and constants
    # The signal level (S) is the sound of the people talking at their location.
    S = 75.0  # dB

    # Step 2: Determine the source strength of the noise sources (train and construction)
    # We assume source strength is the sound level at 1 meter.
    # The formula relating sound level L at distance r is: L(r) = L_1m - 20 * log10(r)

    # For the train: a measurement of 100 dB is taken at 10m from the train and 20m from the people.
    # The sound from the people at that distance (75 - 20*log10(20) ~= 49 dB) is negligible compared to 100 dB.
    # So, we approximate the train's sound level at 10m as 100 dB.
    L_train_10m = 100.0
    # Calculate the train's source strength at 1m.
    L_train_1m = L_train_10m + 20 * math.log10(10)

    # For the construction: a measurement of 115 dB is taken at 20m from construction and 30m from people.
    # The sound from people (75 - 20*log10(30) ~= 45 dB) is negligible compared to 115 dB.
    # So, we approximate the construction's sound level at 20m as 115 dB.
    L_construction_20m = 115.0
    # Calculate the construction's source strength at 1m.
    L_construction_1m = L_construction_20m + 20 * math.log10(20)

    # Step 3: Calculate the noise level from each source at the people's location
    # The train is 30m away, construction is 50m away.
    dist_train = 30.0
    dist_construction = 50.0

    L_train_at_location = L_train_1m - 20 * math.log10(dist_train)
    L_construction_at_location = L_construction_1m - 20 * math.log10(dist_construction)

    # Step 4: Calculate the total noise level (N) at the people's location
    # Intensities must be added, not decibels. I = 10^(L/10)
    intensity_train = 10**(L_train_at_location / 10)
    intensity_construction = 10**(L_construction_at_location / 10)
    total_intensity_noise = intensity_train + intensity_construction

    # Convert total intensity back to decibels. N = 10 * log10(I_total)
    N = 10 * math.log10(total_intensity_noise)

    # Step 5: Calculate the Signal-to-Noise Ratio (SNR)
    SNR = S - N

    # Output the results
    print("--- Calculation Breakdown ---")
    print(f"Signal Level (S) at the location: {S:.2f} dB")
    print(f"Noise from train at the location: {L_train_at_location:.2f} dB")
    print(f"Noise from construction at the location: {L_construction_at_location:.2f} dB")
    print(f"Total Combined Noise Level (N) at the location: {N:.2f} dB")
    print("\n--- Final Equation ---")
    print(f"SNR = Signal - Noise")
    print(f"SNR = {S:.2f} dB - {N:.2f} dB")
    print(f"Final Signal-to-Noise Ratio (SNR): {SNR:.2f} dB")

calculate_snr()
<<<D>>>