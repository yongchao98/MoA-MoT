import cmath
import math

def solve_circuit_voltage():
    """
    Calculates the total DC voltage across the load in the given circuit,
    considering harmonic distortion and parasitic effects.
    """
    # --- Parameters ---
    # Input Signal
    V_in_1 = 1.0  # V, amplitude of the fundamental
    f1 = 915e6    # Hz, fundamental frequency

    # Circuit Components from table and text
    RL = 8e3      # Ohm
    CL = 5e-12    # Farad
    R0 = 50       # Ohm, parasitic resistance base value
    f0 = 915e6    # Hz, reference frequency for parasitic resistance
    C_parasitic = 2e-15 # Farad

    # --- Model Assumptions ---
    # - The circuit acts as an ideal voltage doubler: V_DC = 2 * V_peak_in.
    # - The total DC voltage is the sum of DC voltages from each harmonic.
    # - The system is modeled as a linear voltage divider for each harmonic.

    # --- Calculations ---

    # Total load capacitance is the sum of CL and C_parasitic
    C_total_load = CL + C_parasitic

    harmonics = [1, 3, 5, 7]
    V_in_amplitudes = {}
    V_in_amplitudes[1] = V_in_1
    # Voltage drops by 10% for each higher harmonic
    V_in_amplitudes[3] = V_in_amplitudes[1] * 0.9
    V_in_amplitudes[5] = V_in_amplitudes[3] * 0.9
    V_in_amplitudes[7] = V_in_amplitudes[5] * 0.9

    total_dc_voltage = 0
    dc_voltages_per_harmonic = []

    print("Step-by-step Calculation:")
    print("="*40)
    print(f"Total Load Capacitance (C_L + C_parasitic): {C_total_load*1e12:.4f} pF")
    print("\nProcessing each harmonic component:")

    # Loop through each harmonic
    for n in harmonics:
        f_n = n * f1
        w_n = 2 * math.pi * f_n
        V_in_n = V_in_amplitudes[n]

        # Calculate frequency-dependent series parasitic resistance
        R_s_n = R0 * (f_n / f0)**2

        # Calculate the complex load impedance Z_load = 1 / Y_load
        # Y_load is the admittance of (RL || CL || C_parasitic)
        Y_load_n = 1/RL + 1j * w_n * C_total_load
        Z_load_n = 1 / Y_load_n

        # Use the voltage divider formula to find the voltage at the rectifier input
        V_rect_in_n = V_in_n * (Z_load_n / (R_s_n + Z_load_n))

        # The amplitude is the magnitude of the complex voltage
        V_rect_in_n_mag = abs(V_rect_in_n)

        # Calculate the DC voltage produced by this harmonic using the ideal doubler model
        V_dc_n = 2 * V_rect_in_n_mag
        dc_voltages_per_harmonic.append(V_dc_n)

        # Add to the total DC voltage
        total_dc_voltage += V_dc_n
        
        print(f"\n--- Harmonic n={n} ---")
        print(f"  - Frequency (f_{n}): {f_n/1e6:.2f} MHz")
        print(f"  - Input Voltage Amplitude (V_in_{n}): {V_in_n:.4f} V")
        print(f"  - Series Parasitic Resistance (R_s_{n}): {R_s_n:.4f} Ohm")
        print(f"  - AC Voltage at Rectifier (|V_rect_in_{n}|): {V_rect_in_n_mag:.4f} V")
        print(f"  - Resulting DC Voltage (V_DC_{n}): {V_dc_n:.4f} V")

    print("\n" + "="*40)
    print("Final Calculation Summary:")
    print("The total DC voltage is the sum of the DC contributions from each harmonic:")
    
    equation_symbols = "V_DC_total = " + " + ".join([f"V_DC_{n}" for n in harmonics])
    print(equation_symbols)

    equation_values = "V_DC_total = " + " + ".join([f"{v:.4f}" for v in dc_voltages_per_harmonic])
    print(equation_values)
    
    print(f"\nFinal Total Voltage across C_L = {total_dc_voltage:.4f} V")

# Execute the function
solve_circuit_voltage()