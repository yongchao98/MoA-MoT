import math

def solve_snr():
    """
    Calculates the signal-to-noise ratio for the given sound problem.
    """
    # Signal Level
    L_S = 75  # dB

    # --- Step 1: Calculate the train's noise level at the people's location ---

    # Reference measurement for the train
    L_ref1 = 100  # dB
    d_train_ref1 = 10  # meters
    
    # Target distance for the train
    d_train_final = 30 # meters

    # Calculate train noise level at the people's location (L_N1)
    # L_N1 = L_ref1 - 20 * log10(d_final / d_ref)
    L_N1 = L_ref1 - 20 * math.log10(d_train_final / d_train_ref1)
    
    print("Step 1: Calculate Train Noise at People's Location")
    print(f"Train sound level at 10m: {L_ref1} dB")
    print(f"Train sound level at 30m (L_N1): {L_ref1} - 20*log10({d_train_final}/{d_train_ref1}) = {L_N1:.2f} dB\n")

    # --- Step 2: Calculate the construction's noise level at the people's location ---

    # Reference measurement for the construction
    L_ref2 = 115 # dB
    d_const_ref2 = 20 # meters
    
    # Target distance for the construction
    d_const_final = 50 # meters

    # Calculate construction noise level at the people's location (L_N2)
    L_N2 = L_ref2 - 20 * math.log10(d_const_final / d_const_ref2)

    print("Step 2: Calculate Construction Noise at People's Location")
    print(f"Construction sound level at 20m: {L_ref2} dB")
    print(f"Construction sound level at 50m (L_N2): {L_ref2} - 20*log10({d_const_final}/{d_const_ref2}) = {L_N2:.2f} dB\n")

    # --- Step 3: Combine the noise levels ---
    
    # Convert noise levels (dB) to relative intensities
    I_N1 = 10**(L_N1 / 10)
    I_N2 = 10**(L_N2 / 10)

    # Add intensities
    I_N_total = I_N1 + I_N2

    # Convert total intensity back to dB
    L_N_total = 10 * math.log10(I_N_total)
    
    print("Step 3: Calculate Total Noise Level")
    print(f"Total Noise Level (L_N_total) = 10*log10(10^({L_N1:.2f}/10) + 10^({L_N2:.2f}/10)) = {L_N_total:.2f} dB\n")

    # --- Step 4: Calculate the Signal-to-Noise Ratio (SNR) ---
    
    SNR = L_S - L_N_total
    
    print("Step 4: Calculate Signal-to-Noise Ratio (SNR)")
    print("SNR = Signal Level - Total Noise Level")
    print(f"Final Equation: SNR = {L_S} - {L_N_total:.2f}")
    print(f"Final Answer: {SNR:.2f} dB")

solve_snr()