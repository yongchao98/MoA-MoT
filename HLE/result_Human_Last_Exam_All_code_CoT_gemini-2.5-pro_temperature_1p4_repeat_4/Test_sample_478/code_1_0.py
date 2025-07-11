import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio based on the problem description.
    """
    # Step 1: Define given values
    signal_level = 75.0  # dB, at the people's location

    # Information for calculating noise sources
    train_meas_level = 100.0  # dB
    train_meas_dist = 10.0    # meters
    const_meas_level = 115.0  # dB
    const_meas_dist = 20.0    # meters

    # Distances from noise sources to the people
    train_dist_to_people = 30.0 # meters
    const_dist_to_people = 50.0 # meters

    # Step 2: Calculate the source sound level (L_S at 1m) for the train
    # L(r) = L_S - 20 * log10(r)  =>  L_S = L(r) + 20 * log10(r)
    # As discussed in the plan, we assume the 100 dB is dominated by the train.
    train_source_level = train_meas_level + 20 * math.log10(train_meas_dist)
    
    # Step 3: Calculate the source sound level (L_S at 1m) for the construction
    # Similarly, assume the 115 dB is dominated by the construction.
    const_source_level = const_meas_level + 20 * math.log10(const_meas_dist)

    # Step 4: Calculate the noise level from each source at the people's location
    # L(r) = L_S - 20 * log10(r)
    train_noise_at_location = train_source_level - 20 * math.log10(train_dist_to_people)
    const_noise_at_location = const_source_level - 20 * math.log10(const_dist_to_people)

    # Step 5: Combine the noise levels
    # L_total = 10 * log10(10^(L1/10) + 10^(L2/10))
    total_noise_level = 10 * math.log10(math.pow(10, train_noise_at_location / 10) + math.pow(10, const_noise_at_location / 10))

    # Step 6: Calculate the Signal-to-Noise Ratio (SNR)
    snr = signal_level - total_noise_level

    # Output the final equation and result
    print("Signal level (S) = {:.2f} dB".format(signal_level))
    print("Noise level from train (N_train) = {:.2f} dB".format(train_noise_at_location))
    print("Noise level from construction (N_const) = {:.2f} dB".format(const_noise_at_location))
    print("Total combined noise level (N_total) = {:.2f} dB".format(total_noise_level))
    print("\nFinal SNR Calculation:")
    print("SNR = S - N_total")
    print("{:.2f} - {:.2f} = {:.2f}".format(signal_level, total_noise_level, snr))

solve_snr()
<<<D>>>