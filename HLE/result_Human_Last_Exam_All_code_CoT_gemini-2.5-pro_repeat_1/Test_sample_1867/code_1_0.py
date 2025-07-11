import math

def solve_diode_impedance_matching():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # Step 1: Define Constants and Parameters
    k = 1.380649e-23  # Boltzmann's constant in J/K
    q = 1.602176634e-19 # Elementary charge in C
    T = 300.0           # Ambient temperature in Kelvin
    Io = 1e-9           # Diode reverse saturation current in A
    n = 1.5             # Diode ideality factor
    V1 = 0.78           # Start voltage of linear region in V
    V2 = 0.98           # End voltage of linear region in V
    I2 = 0.445          # Current at V2 in A
    RL_load = 50.0      # Load resistance in ohms
    margin = 0.20       # Startup margin

    # Step 2: Calculate thermal voltage (Vt)
    Vt = (k * T) / q

    # Step 3: Calculate current I1 at V1 using the standard diode equation
    I1 = Io * (math.exp(V1 / (n * Vt)) - 1)

    # Step 4: Calculate the diode's dynamic source resistance (Rs)
    Rs = (V2 - V1) / (I2 - I1)

    # Step 5: Determine the target transformed load resistance (RL_prime) with the startup margin
    # For a negative resistance oscillator to start, the load it sees (RL_prime)
    # must be less than the magnitude of the negative resistance |Rs|.
    Rs_magnitude = abs(Rs)
    RL_prime = Rs_magnitude * (1 - margin)

    # Step 6: Calculate the impedance transformation ratio
    transformation_ratio = RL_prime / RL_load

    # --- Output ---
    print("This script calculates the required impedance transformation ratio for a diode signal source.")
    print("\n--- Intermediate Calculations ---")
    print(f"1. Calculated Current at V1 (I1): {I1:.4f} A")
    print(f"2. Calculated Diode Dynamic Resistance (Rs): {Rs:.4f} Ohms")
    print(f"3. Target Transformed Load Resistance (RL_prime) with {margin*100}% margin: {RL_prime:.4f} Ohms")
    
    print("\n--- Final Equation ---")
    print("The final impedance transformation ratio is calculated as:")
    print(f"Ratio = (Target Transformed Load Resistance) / (Actual Load Resistance)")
    print(f"Ratio = {RL_prime:.4f} / {RL_load:.1f}")
    
    print("\n--- Final Answer ---")
    print(f"The impedance transformation ratio should be: {transformation_ratio:.4f}")
    
    # Return final answer in specified format
    final_answer = f"<<<{transformation_ratio:.4f}>>>"
    print(final_answer)

solve_diode_impedance_matching()