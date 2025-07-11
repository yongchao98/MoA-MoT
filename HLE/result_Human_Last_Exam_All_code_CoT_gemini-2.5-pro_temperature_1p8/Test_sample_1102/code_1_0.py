def solve_welding_parameters():
    """
    Analyzes TIG welding parameters to find the correct current and voltage
    for a turbine blade repair.
    """
    # --- Given Parameters ---
    arc_gap_mm = 6.0
    material = "Inconel 718 Turbine Blade (thin section)"
    travel_speed_mm_s = 0.5
    process = "Manual TIG (GTAW) root pass"

    # --- Step 1: Estimate Arc Voltage ---
    # An empirical model for TIG arc voltage is: V = V_base + (V_gradient * arc_length)
    # Using typical values for Argon shielding gas:
    v_base = 10.0  # Base voltage for anode/cathode drop (V)
    v_gradient = 1.25 # Voltage gradient of the arc plasma (V/mm)

    # Calculate the estimated voltage
    estimated_voltage = v_base + (v_gradient * arc_gap_mm)

    print("Step 1: Calculating the Estimated Arc Voltage")
    print("---------------------------------------------")
    print(f"The arc voltage depends strongly on the arc gap, which is given as {arc_gap_mm} mm.")
    print("Using the physical model: Voltage = Base_Voltage + (Gradient * Arc_Gap)")
    print("Assuming a Base_Voltage of 10 V and a Gradient of 1.25 V/mm:")
    print(f"Voltage = {v_base} V + ({v_gradient} V/mm * {arc_gap_mm} mm)")
    print(f"Calculated Estimated Voltage = {estimated_voltage:.2f} V\n")
    print(f"This result ({estimated_voltage:.2f} V) is very close to 17.5 V, making it the most plausible voltage choice.")

    # --- Step 2: Evaluate Welding Current ---
    print("\nStep 2: Evaluating the Welding Current (Amperage)")
    print("-------------------------------------------------")
    print(f"The application is a precision root pass on a thin '{material}'.")
    print(f"The travel speed is slow ({travel_speed_mm_s} mm/s), which increases heat input per unit length.")
    print("To prevent overheating, distortion, and burn-through on such a delicate component, a low welding current is required.")
    print("A current of 100 A or more would be far too high and destructive.")
    print("A low current like 25 A provides the necessary control for a stable, high-quality repair.\n")

    # --- Conclusion ---
    final_amps = 25
    final_volts = 17.5
    print("\nConclusion: The optimal choice combines the physically correct voltage for the given arc gap with a low, controllable current suitable for the delicate repair.")
    print("\nFinal Recommended Parameters:")
    # The final line prints the full equation with the numbers as requested.
    print(f"Current = {final_amps} A, Voltage = {final_volts} V")


solve_welding_parameters()