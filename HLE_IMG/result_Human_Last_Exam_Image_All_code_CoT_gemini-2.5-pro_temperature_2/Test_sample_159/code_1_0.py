import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Step 1: Define Constants and Initial Values ---
    V_fundamental = 1.0  # V
    f_fundamental = 915e6  # Hz
    R0 = 50.0  # Omega
    f0 = 915e6  # Hz
    RL = 8e3  # Omega

    # --- Step 2: Make Key Assumptions for the Model ---
    # Assume the rectifier's effective input resistance is half the load resistance.
    R_in = RL / 2
    # The voltage drops by 10% for each higher harmonic.
    voltage_drop_factor = 0.9

    harmonics = [1, 3, 5, 7]
    voltages = {}
    parasitic_resistances = {}
    input_powers = {}
    rectifier_powers = {}

    # --- Step 3: Calculate Properties for Each Harmonic ---
    current_v = V_fundamental
    for n in harmonics:
        if n == 1:
            voltages[n] = V_fundamental
        else:
            # Voltage for harmonic n is 90% of the previous one (n-2)
            voltages[n] = voltages[n-2] * voltage_drop_factor

        f_n = n * f_fundamental
        # Parasitic resistance increases with the square of the frequency
        parasitic_resistances[n] = R0 * (f_n / f0)**2

    # --- Step 4: Calculate Power Components for Each Harmonic ---
    total_input_power = 0
    for n in harmonics:
        Vn = voltages[n]
        Rpn = parasitic_resistances[n]
        total_series_resistance = Rpn + R_in
        
        # Power drawn from the source for harmonic n
        # P_in = V_rms^2 / R_total = (V_peak/sqrt(2))^2 / R_total = V_peak^2 / (2 * R_total)
        p_in_n = (Vn**2) / (2 * total_series_resistance)
        input_powers[n] = p_in_n
        total_input_power += p_in_n
        
        # Power delivered to the rectifier for harmonic n
        # P_rect = I_rms^2 * R_in = (V_peak / R_total)^2 / 2 * R_in
        p_rect_n = (Vn**2 * R_in) / (2 * total_series_resistance**2)
        rectifier_powers[n] = p_rect_n

    # --- Step 5: Determine Output Power ---
    # Assume only power from the fundamental is converted to useful DC output
    P_out = rectifier_powers[1]

    # --- Step 6: Calculate Overall Efficiency ---
    efficiency = P_out / total_input_power

    # --- Step 7: Print the Detailed Calculation ---
    print("--- System Efficiency Calculation ---")
    print("This calculation assumes:")
    print(f"1. Rectifier effective input resistance R_in = RL/2 = {R_in:.0f} Ω.")
    print("2. Only power from the fundamental frequency contributes to the DC output.")
    print("-" * 35)

    print("Equation for Efficiency: η = P_out / P_in_total\n")
    
    print(f"P_out = Power delivered to rectifier at fundamental (n=1)")
    V1 = voltages[1]
    Rp1 = parasitic_resistances[1]
    total_R1 = Rp1 + R_in
    print(f"P_out = (V₁² * R_in) / (2 * (Rp₁ + R_in)²) = ({V1:.2f}² * {R_in:.0f}) / (2 * ({Rp1:.0f} + {R_in:.0f})²) = {P_out:.6f} W")
    print("-" * 35)
    
    print(f"P_in_total = Sum of input powers from all harmonics (n=1,3,5,7)")
    print(f"P_in_n = V_n² / (2 * (Rp_n + R_in))")
    for n in harmonics:
        Vn = voltages[n]
        Rpn = parasitic_resistances[n]
        Pin = input_powers[n]
        print(f"P_in_{n} = {Vn:.3f}² / (2 * ({Rpn:.0f} + {R_in:.0f})) = {Pin:.6f} W")
    
    print(f"\nP_in_total = {input_powers[1]:.6f} + {input_powers[3]:.6f} + {input_powers[5]:.6f} + {input_powers[7]:.6f} = {total_input_power:.6f} W")
    print("-" * 35)
    
    print("Final Efficiency Calculation:")
    print(f"η = {P_out:.6f} / {total_input_power:.6f} = {efficiency:.4f}")
    print(f"Overall System Efficiency: {efficiency * 100:.2f}%")
    print(f"\n<<<anwser: {efficiency:.4f}>>>")
    

calculate_efficiency()