import math

def solve_diode_impedance_transformation():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # Step 1: Define all given constants and parameters.
    Io = 1e-9       # Reverse saturation current in Amperes
    n = 1.5         # Diode ideality factor
    T = 300         # Ambient temperature in Kelvin
    V1 = 0.78       # Start voltage of linear region in Volts
    V2 = 0.98       # End voltage of linear region in Volts
    I2 = 0.445      # Current at V2 in Amperes
    R_load = 50.0   # Load resistance in Ohms
    margin = 0.20   # 20% startup margin

    # Physical constants
    k = 1.380649e-23  # Boltzmann constant in J/K
    q = 1.60217663e-19  # Elementary charge in Coulombs

    # Step 2: Calculate the thermal voltage, V_T = kT/q.
    V_T = (k * T) / q
    print(f"Calculated Thermal Voltage (V_T): {V_T:.4f} V")

    # Step 3: Calculate the current I1 at voltage V1 using the diode equation.
    # This is the starting point of the linear signal region.
    try:
        I1 = Io * (math.exp(V1 / (n * V_T)) - 1)
        print(f"Calculated Current at V1 (I1): {I1:.4f} A")
    except OverflowError:
        print("Error: Calculation resulted in overflow. Check input parameters.")
        return

    # Step 4: Calculate the dynamic resistance of the diode, r_d = dV/dI.
    # This acts as the negative source impedance.
    delta_V = V2 - V1
    delta_I = I2 - I1
    if delta_I == 0:
        print("Error: Change in current is zero, cannot calculate dynamic resistance.")
        return
    r_d = delta_V / delta_I
    print(f"Calculated Diode Dynamic Resistance (Source Impedance): {r_d:.4f} Ohms")

    # Step 5: Apply Maximum Power Transfer and the startup margin.
    # The load impedance seen by the diode (R'_L) must be smaller than |r_d| for startup.
    # R'_L = |r_d| * (1 - margin)
    R_L_prime = abs(r_d) * (1 - margin)
    print(f"Required Impedance at Diode with 20% Margin (R'_L): {R_L_prime:.4f} Ohms")

    # Step 6: Determine the impedance transformation ratio.
    # Ratio = Z_load / Z_diode = R_load / R'_L
    if R_L_prime == 0:
        print("Error: Required impedance at the diode is zero, cannot calculate transformation ratio.")
        return
        
    transformation_ratio = R_load / R_L_prime
    print("\nThe impedance transformation ratio is the ratio of the load resistance to the required impedance at the diode.")
    print(f"Final Equation: Ratio = {R_load} Ohms / {R_L_prime:.4f} Ohms")
    print(f"Calculated Impedance Transformation Ratio: {transformation_ratio:.2f}")
    
    # Returning the final answer in the specified format
    final_answer = f"<<<{transformation_ratio:.2f}>>>"
    print(final_answer)

solve_diode_impedance_transformation()