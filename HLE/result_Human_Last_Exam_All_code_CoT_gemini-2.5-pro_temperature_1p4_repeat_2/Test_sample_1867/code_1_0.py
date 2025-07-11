import math

def solve_diode_impedance_ratio():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # --- Given Parameters ---
    Io = 1e-9  # Reverse saturation current in Amps
    n = 1.5    # Diode ideality factor
    T = 300    # Ambient temperature in Kelvin
    V1 = 0.78  # Start voltage of linear region in Volts
    V2 = 0.98  # End voltage of linear region in Volts
    I2 = 0.445 # Current at V2 in Amps
    R_load = 50.0 # Load resistance in Ohms
    margin = 0.20 # Startup margin

    # --- Physical Constants ---
    k = 1.380649e-23  # Boltzmann constant in J/K
    q = 1.60217663e-19  # Elementary charge in Coulombs

    # Step 1: Calculate Thermal Voltage (Vt)
    V_t = (k * T) / q

    # Step 2: Calculate current I1 at voltage V1 using the diode equation
    # This determines the starting point of the linear region.
    # I1 = Io * (e^(V1 / (n * Vt)) - 1)
    exponent = V1 / (n * V_t)
    I1 = Io * (math.exp(exponent) - 1)

    # Step 3: Calculate the dynamic resistance (rd) of the diode
    # rd = dV / dI = (V2 - V1) / (I2 - I1)
    delta_V = V2 - V1
    delta_I = I2 - I1
    # A negative resistance indicates an active source, which is expected.
    r_d = delta_V / delta_I

    # Step 4: The source impedance for matching is the magnitude of rd
    R_source = abs(r_d)

    # Step 5: Apply the 20% startup margin
    R_source_adjusted = R_source * (1 + margin)

    # Step 6: Calculate the required impedance transformation ratio (K)
    # For max power transfer, R_transformed_load = R_source_adjusted
    # K = R_transformed_load / R_load
    K = R_source_adjusted / R_load

    # --- Output the results and the final equation ---
    print(f"Intermediate Calculations:")
    print(f"Thermal Voltage (Vt): {V_t:.4f} V")
    print(f"Calculated Current at V1 (I1): {I1:.4f} A")
    print(f"Dynamic Resistance (rd = dV/dI): {r_d:.4f} Ohms")
    print(f"Source Impedance (|rd|): {R_source:.4f} Ohms")
    print(f"Source Impedance with 20% Margin: {R_source_adjusted:.4f} Ohms")
    print("\nFinal Equation for the Impedance Transformation Ratio (K):")
    print(f"K = (Adjusted Source Impedance) / (Load Resistance)")
    print(f"K = {R_source_adjusted:.4f} Ohms / {R_load:.1f} Ohms")
    print(f"K = {K:.6f}")

solve_diode_impedance_ratio()
<<<0.048663>>>