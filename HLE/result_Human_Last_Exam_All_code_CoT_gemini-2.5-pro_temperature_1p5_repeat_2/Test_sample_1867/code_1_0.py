import math

def solve_impedance_transformation():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """

    # --- Given parameters ---
    Io = 1e-9  # Reverse saturation current in Amperes
    n = 1.5    # Diode ideality factor
    T = 300    # Ambient temperature in Kelvin
    V1 = 0.78  # Start voltage of linear region in Volts
    V2 = 0.98  # End voltage of linear region in Volts
    I2 = 0.445 # Current at V2 in Amperes
    R_load = 50.0  # Load resistance in Ohms
    margin = 0.20  # 20% startup margin

    # --- Physical constants ---
    k = 1.380649e-23  # Boltzmann constant in J/K
    q = 1.60217663e-19 # Elementary charge in Coulombs

    # Step 1: Calculate the thermal voltage (VT)
    V_T = (k * T) / q
    print(f"Step 1: Calculated Thermal Voltage (VT)")
    print(f"VT = (k * T) / q = ({k:.4e} * {T}) / {q:.4e} = {V_T:.4f} V\n")

    # Step 2: Calculate the current I1 at the start of the linear region
    # using the standard diode equation.
    I1 = Io * (math.exp(V1 / (n * V_T)) - 1)
    print(f"Step 2: Calculated Current I1 at V1 = {V1} V")
    print(f"I1 = Io * (exp(V1 / (n * VT)) - 1)")
    print(f"I1 = {Io:.1e} * (exp({V1} / ({n} * {V_T:.4f})) - 1) = {I1:.4f} A\n")

    # Step 3: Calculate the dynamic source resistance (R_source) of the diode
    # in the specified linear region.
    delta_V = V2 - V1
    delta_I = I2 - I1
    R_source = delta_V / delta_I
    print(f"Step 3: Calculated Dynamic Source Resistance (R_source)")
    print(f"R_source = (V2 - V1) / (I2 - I1)")
    print(f"R_source = ({V2} - {V1}) / ({I2} - {I1:.4f}) = {R_source:.4f} Ohms\n")

    # Step 4: Determine the target transformed load resistance with the startup margin
    # The magnitude of the source resistance must be >= (1+margin) * transformed_load
    R_transformed_load = abs(R_source) / (1 + margin)
    print(f"Step 4: Calculated Target Transformed Load Resistance (R_transformed_load)")
    print(f"R_transformed_load = |R_source| / (1 + margin)")
    print(f"R_transformed_load = {abs(R_source):.4f} / (1 + {margin}) = {R_transformed_load:.4f} Ohms\n")

    # Step 5: Calculate the final impedance transformation ratio
    transformation_ratio = R_load / R_transformed_load
    print(f"Step 5: Final Impedance Transformation Ratio Calculation")
    print(f"Ratio = R_load / R_transformed_load")
    print(f"Ratio = {R_load} / {R_transformed_load:.4f} = {transformation_ratio:.4f}\n")
    
    print("Final Answer:")
    print(f"The impedance transformation ratio should be {transformation_ratio:.2f}:1")
    
    # Final answer in the required format
    print(f"\n<<<{transformation_ratio:.2f}>>>")

solve_impedance_transformation()