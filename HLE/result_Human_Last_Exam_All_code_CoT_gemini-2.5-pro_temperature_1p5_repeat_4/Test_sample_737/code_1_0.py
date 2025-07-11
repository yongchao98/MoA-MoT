import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Parameters ---
    # Convert all units to SI (meters, Newtons, Pascals, kg) for consistency.

    # Plate Bearing Test Details
    P_plate = 30 * 1000  # Load in Newtons (30 kN)
    d_plate = 305 / 1000 # Plate diameter in meters (305 mm)
    a_plate = d_plate / 2 # Plate radius in meters

    # Subgrade Test Results (Layer 2)
    delta_subgrade = 2460 / 1e6 # Deflection in meters (2460 μm)

    # Trial Pavement Test Results (Layer 1 on Layer 2)
    h_trial = 300 / 1000 # Trial pavement thickness in meters (300 mm)
    delta_trial = 1080 / 1e6 # Deflection in meters (1080 μm)

    # Design Parameters and Limits
    max_deflection_design = 1.00 / 1000 # Max allowable deflection in meters (1.00 mm)
    design_load_mass = 1.80 * 1000 # Design wheel load in kg (1.80 ton)
    g = 9.81 # Acceleration due to gravity in m/s^2
    P_design = design_load_mass * g # Design wheel load in Newtons
    p_design = 600 * 1000 # Tyre pressure in Pascals (600 kN/m^2)

    # --- Step-by-step Solution ---

    print("Step 1: Determine the properties of the subgrade (Layer 2).")

    # Pressure under the plate during the test
    p_plate = P_plate / (math.pi * a_plate**2)
    print(f"The pressure from the plate bearing test (p_plate) is {p_plate:.2f} Pa.")

    # Boussinesq's equation for deflection on a single layer is used to find E2.
    # For a flexible plate with mu=0.5, Δ = 1.5 * p * a / E.
    # Rearranging for E2: E2 = 1.5 * p_plate * a_plate / delta_subgrade
    E2 = (1.5 * p_plate * a_plate) / delta_subgrade
    print(f"The elastic modulus of the subgrade (E2) is calculated to be {E2 / 1e6:.2f} MPa.\n")

    print("Step 2: Determine the parameters for the design wheel load.")
    # The radius of the design wheel load contact area 'a_design' is found from the pressure and load.
    # p_design = P_design / (pi * a_design^2) => a_design = sqrt(P_design / (pi * p_design))
    a_design = math.sqrt(P_design / (math.pi * p_design))
    print(f"The radius of the design wheel load (a_design) is {a_design:.4f} m.\n")

    print("Step 3: Calculate and compare the deflection factor F2 for both cases.")
    # For a two-layer system, deflection Δ = (1.5 * p * a / E2) * F2.
    # We can calculate the required F2 factor for the design case.
    F2_design = (max_deflection_design * E2) / (1.5 * p_design * a_design)
    print(f"The required deflection factor F2 for the design (F2_design) is {F2_design:.4f}.")

    # And the F2 factor from the trial section case.
    F2_trial = (delta_trial * E2) / (1.5 * p_plate * a_plate)
    print(f"The deflection factor F2 from the trial section (F2_trial) is {F2_trial:.4f}.")
    print("Observation: The F2 values are almost identical.\n")

    print("Step 4: Calculate the final required thickness.")
    print("Since F2 is a function of h/a and E1/E2, and the F2 values are equal, we can assume the h/a ratios are also equal.")
    print("Therefore: h_design / a_design = h_trial / a_plate")
    print("Rearranging for h_design gives:")

    # Solve for the required pavement thickness, h_design
    h_design = h_trial * (a_design / a_plate)

    print(f"h_design = {h_trial:.3f} m * ({a_design:.4f} m / {a_plate:.4f} m)")

    # Convert the final result to mm.
    h_design_mm = h_design * 1000

    print(f"\nThe required pavement thickness (h_design) is {h_design_mm:.2f} mm.")

    # --- Final Answer ---
    print(f"<<<{h_design_mm:.2f}>>>")

solve_pavement_thickness()