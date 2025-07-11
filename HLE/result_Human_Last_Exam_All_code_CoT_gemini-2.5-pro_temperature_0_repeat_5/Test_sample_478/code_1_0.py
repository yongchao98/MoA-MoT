import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio for the given scenario.
    """
    # Step 1: Define signal level and other given values.
    # The signal is the sound from the people at their location.
    L_signal = 75.0  # dB

    # Measurement 1: Train and People
    L_meas_1 = 100.0  # dB
    r_train_meas_1 = 10.0  # meters

    # Measurement 2: Construction and People
    L_meas_2 = 115.0  # dB
    r_const_meas_2 = 20.0  # meters

    # Distances of noise sources from the final location (where the people are)
    r_train_final = 30.0  # meters
    r_const_final = 50.0  # meters

    print("This script calculates the Signal-to-Noise Ratio (SNR) at the people's location.")
    print("-" * 50)
    print(f"Step 1: The signal level from the people is given as {L_signal} dB.")
    print("-" * 50)

    # Step 2: Determine the effective noise levels at the measurement points.
    # We assume the sound from the people is negligible at the measurement points,
    # so the measured level is dominated by the primary noise source.
    L_train_at_meas_1 = L_meas_1
    L_const_at_meas_2 = L_meas_2
    print("Step 2: Determine the sound level of each noise source from the measurements.")
    print(f"The sound level of the train at {r_train_meas_1}m is taken as {L_train_at_meas_1} dB.")
    print(f"The sound level of the construction at {r_const_meas_2}m is taken as {L_const_at_meas_2} dB.")
    print("-" * 50)

    # Step 3: Calculate the sound level of each noise source at the people's location.
    # The formula for sound level change with distance is: L2 = L1 - 20 * log10(r2 / r1)
    L_train_final = L_train_at_meas_1 - 20 * math.log10(r_train_final / r_train_meas_1)
    L_const_final = L_const_at_meas_2 - 20 * math.log10(r_const_final / r_const_meas_2)
    print("Step 3: Calculate the noise level from each source at the people's location.")
    print(f"Train's noise level at {r_train_final}m = {L_train_at_meas_1:.2f} - 20*log10({r_train_final}/{r_train_meas_1}) = {L_train_final:.2f} dB")
    print(f"Construction's noise level at {r_const_final}m = {L_const_at_meas_2:.2f} - 20*log10({r_const_final}/{r_const_meas_2}) = {L_const_final:.2f} dB")
    print("-" * 50)

    # Step 4: Combine the two noise levels to find the total noise level.
    # Convert dB to relative intensity (I_rel = 10^(L/10)), add intensities, and convert back to dB.
    I_rel_train = 10**(L_train_final / 10)
    I_rel_const = 10**(L_const_final / 10)
    I_rel_total_noise = I_rel_train + I_rel_const
    L_noise_total = 10 * math.log10(I_rel_total_noise)
    print("Step 4: Combine the noise levels by adding their intensities.")
    print(f"Total Noise Level = 10*log10(10^({L_train_final:.2f}/10) + 10^({L_const_final:.2f}/10)) = {L_noise_total:.2f} dB")
    print("-" * 50)

    # Step 5: Calculate the final Signal-to-Noise Ratio (SNR).
    snr = L_signal - L_noise_total
    print("Step 5: Calculate the final Signal-to-Noise Ratio (SNR).")
    print("SNR = Signal Level - Total Noise Level")
    print(f"SNR = {L_signal:.2f} dB - {L_noise_total:.2f} dB = {snr:.2f} dB")
    print("-" * 50)

solve_snr()