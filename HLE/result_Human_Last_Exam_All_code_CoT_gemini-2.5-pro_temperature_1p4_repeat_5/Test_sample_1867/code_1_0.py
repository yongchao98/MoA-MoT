import math

def solve_diode_impedance():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # Given parameters
    Io = 1e-9      # Reverse saturation current (A)
    n = 1.5        # Diode ideality factor
    T = 300        # Ambient temperature (K)
    V1 = 0.78      # Start voltage of linear region (V)
    V2 = 0.98      # End voltage of linear region (V)
    I2 = 0.445     # Current at V2 (A)
    R_load = 50.0  # Load resistance (ohms)
    margin = 0.20  # Startup margin

    # Physical constants
    k = 1.380649e-23  # Boltzmann constant (J/K)
    q = 1.602176e-19  # Elementary charge (C)

    # --- Step 1: Calculate Diode Current at V1 ---
    Vt = (k * T) / q  # Thermal voltage
    I1 = Io * (math.exp(V1 / (n * Vt)) - 1)
    
    print(f"Step 1: Calculated current I1 at V1 = {V1:.2f} V is {I1:.4f} A")

    # --- Step 2: Determine Dynamic Resistance ---
    delta_V = V2 - V1
    delta_I = I2 - I1
    Rd = delta_V / delta_I
    
    print(f"Step 2: The diode's dynamic resistance (source impedance) is:")
    print(f"        Rd = ({V2:.2f} - {V1:.2f}) V / ({I2:.4f} - {I1:.4f}) A = {Rd:.4f} ohms")
    
    # The magnitude of the resistance is the impedance we need to match
    Rd_mag = abs(Rd)
    print(f"        The magnitude of the source impedance |Rd| is {Rd_mag:.4f} ohms")

    # --- Step 3: Apply Startup Margin ---
    # The target impedance seen by the diode must be less than |Rd| for startup
    R_target = Rd_mag * (1 - margin)
    
    print(f"\nStep 3: Applying a {margin*100}% startup margin, the target impedance is:")
    print(f"        R_target = {Rd_mag:.4f} * (1 - {margin:.2f}) = {R_target:.4f} ohms")

    # --- Step 4: Calculate Impedance Transformation Ratio ---
    # The ratio is the target impedance divided by the actual load impedance
    transformation_ratio = R_target / R_load
    
    print(f"\nStep 4: The impedance transformation ratio from load to diode is:")
    print(f"        Ratio = R_target / R_load")
    print(f"        Ratio = {R_target:.4f} / {R_load:.1f} = {transformation_ratio:.6f}")
    
    # Final answer
    print(f"\nThe final impedance transformation ratio should be {transformation_ratio:.6f}")
    return transformation_ratio

# Execute the function and print the final answer in the required format.
final_answer = solve_diode_impedance()
print(f"<<<{final_answer:.6f}>>>")
