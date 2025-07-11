import math

def solve_impedance_transformation():
    """
    Calculates the optimal impedance transformation ratio for a diode signal source.
    """
    # --- Given parameters ---
    Io = 1e-9      # Reverse saturation current in Amperes
    n_ideality = 1.5 # Diode ideality factor
    T = 300        # Ambient temperature in Kelvin
    V1 = 0.78      # Start voltage of the linear region in Volts
    V2 = 0.98      # End voltage of the linear region in Volts
    I2 = 0.445     # Current at V2 in Amperes
    RL = 50.0      # Load resistance in Ohms
    margin = 0.20  # Startup margin

    # --- Physical constants ---
    k = 1.380649e-23  # Boltzmann's constant in J/K
    q = 1.60217663e-19  # Elementary charge in Coulombs

    # Step 1: Calculate the thermal voltage (Vt)
    Vt = (k * T) / q
    print(f"Step 1: The thermal voltage (Vt) at {T} K is {Vt:.5f} V.")

    # Step 2: Calculate the current I1 at V1 using the Shockley diode equation
    # I = Io * (exp(V / (n * Vt)) - 1)
    # The '-1' is negligible but included for completeness.
    I1 = Io * (math.exp(V1 / (n_ideality * Vt)) - 1)
    print(f"Step 2: The current (I1) at the start of the linear region (V1 = {V1} V) is {I1:.5f} A.")

    # Step 3: Calculate the dynamic resistance (rd) of the diode in the linear region
    # This is the source impedance (Rs).
    delta_V = V2 - V1
    delta_I = I2 - I1
    rd = delta_V / delta_I
    print(f"Step 3: The dynamic resistance (rd) of the diode source is calculated as dV/dI.")
    print(f"       rd = ({V2} V - {V1} V) / ({I2} A - {I1:.5f} A) = {rd:.5f} Ohms.")
    print("       (Note: The negative resistance indicates the diode is acting as a power source).")

    # Step 4: Determine the ideal impedance transformation ratio (K_ideal)
    # For max power transfer from a negative resistance source (-Rs) to a load (RL'), RL' = Rs.
    # The transformed load impedance is K * RL. So, K * RL = -rd.
    K_ideal = -rd / RL
    print(f"\nStep 4: For maximum power transfer, the transformed load must equal -rd.")
    print(f"       Ideal ratio K = -rd / RL = {-rd:.5f} Ohms / {RL} Ohms = {K_ideal:.5f}.")

    # Step 5: Apply the 20% startup margin
    K_final = K_ideal * (1 + margin)
    print(f"\nStep 5: Applying the {margin*100}% startup margin.")
    print(f"       Final Ratio = {K_ideal:.5f} * (1 + {margin}) = {K_final:.5f}.")

    # --- Final Output ---
    print("\n--- Final Equation ---")
    print("The impedance transformation ratio is calculated as:")
    print(f"Ratio = ((- (V2 - V1) / (I2 - I1)) / RL) * (1 + Margin)")
    # Print the equation with all the intermediate values
    final_equation_str = (
        f"Ratio = ((- ({V2} - {V1}) / ({I2} - {I1:.5f})) / {RL}) * (1 + {margin})"
    )
    print(final_equation_str)
    print(f"Ratio = {K_final:.5f}")
    
    # Final answer in the required format
    print(f"\n<<<{K_final:.5f}>>>")


solve_impedance_transformation()