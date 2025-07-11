import math

def calculate_impedance_ratio():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # Given parameters
    Io = 10**-9  # Reverse saturation current in A
    n = 1.5      # Diode ideality factor
    T = 300      # Ambient temperature in K
    V1 = 0.78    # Start voltage of linear region in V
    V2 = 0.98    # End voltage of linear region in V
    I2 = 0.445   # Current at V2 in A
    RL = 50      # Load resistance in ohms
    margin = 0.20 # Startup margin (20%)

    # Physical constants
    k = 1.380649e-23 # Boltzmann constant in J/K
    q = 1.60217663e-19 # Elementary charge in C

    # --- Plan Explanation ---
    print("This script calculates the impedance transformation ratio based on the following plan:")
    print("1. Calculate the thermal voltage (Vt) using fundamental constants and temperature.")
    print("2. Calculate the diode current (I1) at V1 using the Shockley diode equation, which defines the start of the linear operating region.")
    print("3. Calculate the diode's dynamic resistance (rd) from the slope of the linear region between (V1, I1) and (V2, I2).")
    print("4. Determine the required transformed load resistance (RL_transformed) that satisfies the 20% startup margin requirement.")
    print("5. Calculate the final impedance transformation ratio by dividing the load resistance by the required transformed resistance.")
    print("-" * 30)

    # --- Calculations ---

    # Step 1: Calculate thermal voltage Vt
    Vt = (k * T) / q
    print(f"Step 1: Calculate Thermal Voltage (Vt)")
    print(f"Vt = (k * T) / q = ({k:.4e} J/K * {T} K) / {q:.4e} C = {Vt:.4f} V")
    print("-" * 30)

    # Step 2: Calculate current I1 at V1
    # The term "-1" is negligible at this forward bias but included for completeness.
    I1 = Io * (math.exp(V1 / (n * Vt)) - 1)
    print(f"Step 2: Calculate Current at V1 (I1)")
    print(f"I1 = Io * (exp(V1 / (n * Vt)) - 1)")
    print(f"I1 = {Io:.1e} A * (exp({V1} V / ({n} * {Vt:.4f} V)) - 1) = {I1:.4f} A")
    print("-" * 30)

    # Step 3: Calculate dynamic resistance rd
    rd = (V2 - V1) / (I2 - I1)
    print(f"Step 3: Calculate Dynamic Resistance (rd)")
    print(f"rd = (V2 - V1) / (I2 - I1)")
    print(f"rd = ({V2} V - {V1} V) / ({I2} A - {I1:.4f} A) = {rd:.4f} Ohms")
    print("Note: The negative resistance indicates the diode is an active signal source.")
    print("-" * 30)
    
    # Step 4: Calculate the required transformed load resistance
    RL_transformed = abs(rd) / (1 + margin)
    print(f"Step 4: Calculate Required Transformed Load Resistance (RL_transformed)")
    print(f"RL_transformed = |rd| / (1 + margin)")
    print(f"RL_transformed = |{rd:.4f}| Ohms / (1 + {margin}) = {RL_transformed:.4f} Ohms")
    print("-" * 30)

    # Step 5: Calculate the impedance transformation ratio
    ratio = RL / RL_transformed
    print(f"Step 5: Calculate Impedance Transformation Ratio")
    print("The ratio is defined as Load Impedance / Required Transformed Load Impedance.")
    print("\nFinal Equation:")
    print(f"Ratio = {RL} / ({abs(rd):.4f} / (1 + {margin})) = {ratio:.2f}")
    
    # --- Final Answer ---
    print(f"\n<<< {ratio:.2f} >>>")

calculate_impedance_ratio()