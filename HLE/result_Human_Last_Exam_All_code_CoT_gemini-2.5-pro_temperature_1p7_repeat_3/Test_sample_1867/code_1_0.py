import math

def calculate_impedance_transformation_ratio():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # --- Given Parameters ---
    Io = 1e-9  # Reverse saturation current in Amperes
    n = 1.5    # Diode ideality factor
    T = 300    # Ambient temperature in Kelvin
    V1 = 0.78  # Start voltage in Volts
    V2 = 0.98  # End voltage in Volts
    I2 = 0.445 # Current at V2 in Amperes
    RL_load = 50 # Load resistance in Ohms
    margin = 0.20 # Startup margin (20%)

    # --- Physical Constants ---
    k = 1.380649e-23  # Boltzmann constant in J/K
    q = 1.60217663e-19  # Electron charge in Coulombs

    # --- Step 1: Calculate Thermal Voltage (Vt) ---
    Vt = (k * T) / q

    # --- Step 2: Calculate Current I1 at V1 using the diode equation ---
    # I1 = Io * (exp(V1 / (n * Vt)) - 1)
    # The '-1' term is negligible for this forward bias voltage but included for completeness.
    exponent = V1 / (n * Vt)
    I1 = Io * (math.exp(exponent) - 1)

    # --- Step 3: Calculate the Diode's Dynamic Resistance (Rs) ---
    delta_V = V2 - V1
    delta_I = I2 - I1
    Rs = delta_V / delta_I
    Rs_mag = abs(Rs)

    # --- Step 4: Calculate the required transformed impedance for startup ---
    # The transformer must present this impedance to the diode.
    # RL_transformed = |Rs| / (1 + margin)
    RL_transformed = Rs_mag / (1 + margin)

    # --- Step 5: Calculate the impedance transformation ratio ---
    # The ratio is defined as the ratio of the load impedance to the
    # impedance seen by the diode.
    transformation_ratio = RL_load / RL_transformed

    # --- Output the results step-by-step ---
    print(f"Step 1: The diode's dynamic resistance (source impedance) is being calculated.")
    print(f"         Voltage changes from V1 = {V1:.2f} V to V2 = {V2:.2f} V.")
    print(f"         Current at V1 is calculated to be I1 = {I1:.4f} A.")
    print(f"         Given current at V2 is I2 = {I2:.4f} A.")
    print(f"         Dynamic Resistance Rs = dV/dI = ({V2} - {V1}) / ({I2:.4f} - {I1:.4f}) = {Rs:.4f} Ohms.")

    print(f"\nStep 2: The target impedance for the diode is calculated with a {margin*100:.0f}% startup margin.")
    print(f"         Target Impedance = |Rs| / (1 + margin) = {Rs_mag:.4f} / (1 + {margin:.2f}) = {RL_transformed:.4f} Ohms.")

    print(f"\nStep 3: The impedance transformation ratio is calculated.")
    print(f"         The transformer must match the {RL_load:.1f} Ohm load to the target impedance of {RL_transformed:.4f} Ohms.")
    print(f"\nFinal Equation:")
    print(f"Transformation Ratio = Load Impedance / Diode-Side Impedance")
    print(f"Transformation Ratio = {RL_load:.1f} / {RL_transformed:.4f} = {transformation_ratio:.4f}")

    return transformation_ratio

# Execute the function and get the final answer
final_answer = calculate_impedance_transformation_ratio()
print(f"\n<<<The final impedance transformation ratio is {final_answer:.4f}>>>")
