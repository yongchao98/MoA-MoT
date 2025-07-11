import math

def calculate_transformation_ratio():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # --- Step 1: Define constants and parameters ---
    k = 1.380649e-23  # Boltzmann's constant in J/K
    q = 1.60217663e-19  # Elementary charge in C
    T = 300.0          # Ambient temperature in K
    Io = 1.0e-9        # Reverse saturation current in A
    n = 1.5            # Diode ideality factor
    V1 = 0.78          # Start voltage of linear region in V
    V2 = 0.98          # End voltage of linear region in V
    I2 = 0.445         # Current at V2 in A
    R_load = 50.0      # Load resistance in ohms
    margin = 0.20      # Startup margin (20%)

    # --- Step 2: Calculate thermal voltage (Vt) ---
    Vt = (k * T) / q
    
    # --- Step 3: Calculate current I1 using the Shockley diode equation ---
    exponent = V1 / (n * Vt)
    I1 = Io * (math.exp(exponent) - 1)
    
    # --- Step 4: Calculate the dynamic source resistance (Rs) ---
    # Rs will be negative, indicating an active device
    delta_V = V2 - V1
    delta_I = I2 - I1
    Rs = delta_V / delta_I
    
    # --- Step 5: Calculate target impedance with startup margin ---
    # For oscillator startup, R_target > |Rs|.
    R_target = abs(Rs) * (1 + margin)
    
    # --- Step 6: Calculate the impedance transformation ratio ---
    transformation_ratio = R_target / R_load
    
    print("This script calculates the required impedance transformation ratio for a diode source.")
    print("-" * 50)
    print(f"Calculated dynamic source resistance (Rs): {Rs:.4f} Ohms")
    print(f"Magnitude of source resistance |Rs|: {abs(Rs):.4f} Ohms")
    print(f"Target impedance for diode with {margin*100}% startup margin (|Rs| * {1+margin}): {R_target:.4f} Ohms")
    print(f"Final Load Resistance (R_load): {R_load:.1f} Ohms")
    print("-" * 50)
    print("The impedance transformation ratio is the target impedance divided by the load resistance.")
    print(f"Transformation Ratio = {R_target:.4f} / {R_load:.1f} = {transformation_ratio:.4f}")
    
    return transformation_ratio

# Execute the calculation and store the final answer
final_answer = calculate_transformation_ratio()
print(f"\n<<<final_answer>>>\n{final_answer:.5f}\n<<<final_answer>>>")