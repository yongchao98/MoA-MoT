import math

def solve_snr():
    """
    This function calculates the signal-to-noise ratio based on the problem description.
    """
    # --- Step 1: Define Signal and Noise Source Information ---
    signal_db = 75.0

    # Train noise reference
    L_train_ref = 100.0
    r_train_ref = 10.0
    r_train_to_people = 30.0

    # Construction noise reference
    L_const_ref = 115.0
    r_const_ref = 20.0
    r_const_to_people = 50.0

    print("Step-by-step Calculation:\n")
    print(f"1. The signal level at the people's location is S = {signal_db} dB.\n")

    # --- Step 2: Calculate Noise Level of Each Source at the People's Location ---
    # Formula: L2 = L1 - 20 * log10(r2 / r1)

    # Train
    L_train_at_people = L_train_ref - 20 * math.log10(r_train_to_people / r_train_ref)
    print(f"2. Calculate noise from the train at the people's location (30m):")
    print(f"   L_train = {L_train_ref:.2f} dB - 20 * log10({r_train_to_people:.0f}m / {r_train_ref:.0f}m) = {L_train_at_people:.2f} dB\n")

    # Construction
    L_const_at_people = L_const_ref - 20 * math.log10(r_const_to_people / r_const_ref)
    print(f"3. Calculate noise from construction at the people's location (50m):")
    print(f"   L_construction = {L_const_ref:.2f} dB - 20 * log10({r_const_to_people:.0f}m / {r_const_ref:.0f}m) = {L_const_at_people:.2f} dB\n")

    # --- Step 3: Combine Noise Sources ---
    # Convert dB to relative intensity: I_rel = 10^(L/10)
    I_train_rel = 10**(L_train_at_people / 10)
    I_const_rel = 10**(L_const_at_people / 10)

    # Add intensities
    I_total_noise_rel = I_train_rel + I_const_rel

    # Convert total intensity back to dB: L_total = 10 * log10(I_total)
    L_total_noise = 10 * math.log10(I_total_noise_rel)
    print("4. Combine the two noise sources:")
    print(f"   - Convert each noise level to relative intensity.")
    print(f"   - Add the intensities together.")
    print(f"   - Convert the total intensity back to decibels.")
    print(f"   Total Noise Level (N) = {L_total_noise:.2f} dB\n")

    # --- Step 4: Calculate Signal-to-Noise Ratio (SNR) ---
    snr = signal_db - L_total_noise
    print("5. Calculate the Signal-to-Noise Ratio (SNR):")
    print("   SNR = Signal (dB) - Total Noise (dB)")
    print(f"   SNR = {signal_db:.2f} dB - {L_total_noise:.2f} dB\n")
    print(f"Final SNR = {snr:.2f} dB")


if __name__ == '__main__':
    solve_snr()