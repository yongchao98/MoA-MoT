import math

def calculate_snr():
    """
    Calculates the signal-to-noise ratio at the location of the people.
    """
    # --- Given Information ---
    # Signal level (S) at the people's location
    S = 75.0  # dB

    # Distances of noise sources from the people
    dist_train_to_people = 30.0  # meters
    dist_const_to_people = 50.0  # meters

    # Measurement point 1: Train and People
    L_train_at_10m = 100.0  # dB (Approximating, as people's noise is negligible here)
    dist_train_ref = 10.0   # meters

    # Measurement point 2: Construction and People
    L_const_at_20m = 115.0  # dB (Approximating, as people's noise is negligible here)
    dist_const_ref = 20.0   # meters

    # --- Step 1: Calculate noise level from each source at the people's location ---
    # The formula for sound level change with distance is: L2 = L1 - 20 * log10(r2 / r1)

    # Noise from train at the people's location (r2 = 30m)
    N_train_at_people = L_train_at_10m - 20 * math.log10(dist_train_to_people / dist_train_ref)

    # Noise from construction at the people's location (r2 = 50m)
    N_const_at_people = L_const_at_20m - 20 * math.log10(dist_const_to_people / dist_const_ref)

    # --- Step 2: Calculate the total noise level ---
    # Convert dB to intensity-like values, sum them, and convert back to dB.
    # I_rel = 10^(L/10)
    I_train_rel = 10**(N_train_at_people / 10)
    I_const_rel = 10**(N_const_at_people / 10)
    
    # Total relative intensity of the noise
    I_total_noise_rel = I_train_rel + I_const_rel
    
    # Convert total intensity back to dB
    N_total = 10 * math.log10(I_total_noise_rel)

    # --- Step 3: Calculate the Signal-to-Noise Ratio (SNR) ---
    SNR = S - N_total
    
    # --- Final Output ---
    print(f"The signal level at the location is {S:.2f} dB.")
    print(f"The noise level from the train at the location is {N_train_at_people:.2f} dB.")
    print(f"The noise level from the construction at the location is {N_const_at_people:.2f} dB.")
    print(f"The total combined noise level at the location is {N_total:.2f} dB.")
    print("\n--- Final Calculation ---")
    print(f"SNR (dB) = Signal (dB) - Total Noise (dB)")
    print(f"SNR = {S:.2f} - {N_total:.2f} = {SNR:.2f} dB")


if __name__ == "__main__":
    calculate_snr()
<<<D>>>