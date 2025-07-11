import math

def calculate_snr():
    """
    Calculates the signal-to-noise ratio for the given problem.
    """
    # Step 1: Define the signal level and known noise source parameters
    L_signal = 75.0  # dB, signal at the location

    L_train_ref = 100.0  # dB, train noise at reference distance
    r_train_ref = 10.0   # meters, reference distance for train
    r_train_target = 30.0 # meters, distance from train to people

    L_construction_ref = 115.0 # dB, construction noise at reference distance
    r_construction_ref = 20.0  # meters, reference distance for construction
    r_construction_target = 50.0 # meters, distance from construction to people

    print("Step 1: Identify the signal level.")
    print(f"Signal Level (L_Signal) = {L_signal:.2f} dB\n")

    # Step 2: Calculate the noise level from each source at the people's location
    print("Step 2: Calculate the noise level from each source at the target location.")
    # For the train
    L_train_target = L_train_ref - 20 * math.log10(r_train_target / r_train_ref)
    print(f"Noise from train at {r_train_target:.0f}m = {L_train_ref:.2f} dB - 20*log10({r_train_target:.0f}/{r_train_ref:.0f}) = {L_train_target:.2f} dB")

    # For the construction
    L_construction_target = L_construction_ref - 20 * math.log10(r_construction_target / r_construction_ref)
    print(f"Noise from construction at {r_construction_target:.0f}m = {L_construction_ref:.2f} dB - 20*log10({r_construction_target:.0f}/{r_construction_ref:.0f}) = {L_construction_target:.2f} dB\n")

    # Step 3: Combine the noise levels
    print("Step 3: Combine the noise levels by converting to intensity, adding, and converting back to dB.")
    # Convert dB to relative intensity (I = 10^(L/10))
    I_train_target = 10**(L_train_target / 10)
    I_construction_target = 10**(L_construction_target / 10)

    # Add intensities
    I_total_noise = I_train_target + I_construction_target

    # Convert total intensity back to dB (L = 10*log10(I))
    L_total_noise = 10 * math.log10(I_total_noise)
    print(f"Total Noise Level (L_Noise) = 10*log10(10^({L_train_target:.2f}/10) + 10^({L_construction_target:.2f}/10)) = {L_total_noise:.2f} dB\n")

    # Step 4: Calculate the Signal-to-Noise Ratio (SNR)
    print("Step 4: Calculate the final Signal-to-Noise Ratio (SNR).")
    snr = L_signal - L_total_noise
    print(f"SNR = L_Signal - L_Noise")
    print(f"SNR = {L_signal:.2f} dB - {L_total_noise:.2f} dB = {snr:.2f} dB")

calculate_snr()
<<<D>>>