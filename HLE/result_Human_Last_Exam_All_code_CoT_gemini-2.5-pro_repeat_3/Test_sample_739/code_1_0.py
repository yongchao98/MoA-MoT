import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a road surface by iteratively
    finding the design rainfall intensity and then applying Manning's equation.
    """
    # --- Step 1: Define Road and Pavement Parameters ---
    # Number of lanes per direction
    num_lanes = 3
    # Lane width in meters
    lane_width_m = 3.6
    # Cross-slope (S) as a decimal
    cross_slope_S = 1.75 / 100
    # Manning's roughness coefficient (n) for rough-textured asphalt.
    # This is an assumed standard value as it was not provided.
    manning_n = 0.016

    # The total width of the traveled way is the drainage path length (L)
    drainage_path_L_m = num_lanes * lane_width_m

    # --- Step 2: Determine Design Rainfall Intensity (i) via Iteration ---

    # Since Rainfall-Duration-Frequency (IDF) curves were not provided, we assume
    # a standard IDF relationship for a 10-year return period storm.
    # The formula is i = a / (t_d + b)^c, with parameters converted to SI units
    # (i in mm/hr, t_d in minutes).
    idf_a = 2026.9
    idf_b = 11.9
    idf_c = 0.89

    # The time of concentration (t_c) depends on rainfall intensity (i), and
    # the design intensity (i) depends on storm duration (t_c). We iterate
    # to find a stable solution.

    def calculate_tc(intensity_i, n, L, S):
        """Calculates time of concentration (t_c) in minutes using the FHWA HEC-22 formula."""
        K_si = 3.28
        numerator = K_si * (n * L)**0.6
        denominator = (intensity_i**0.4) * (S**0.3)
        return numerator / denominator

    def calculate_i(duration_t):
        """Calculates rainfall intensity (i) in mm/hr from the IDF curve."""
        return idf_a / ((duration_t + idf_b)**idf_c)

    # Perform the iterative calculation
    design_i_mmhr = 150.0  # Initial guess for intensity in mm/hr
    tolerance = 0.1
    print("This script calculates the design water film thickness on a road.")
    print("It assumes a 10-year storm event and standard pavement roughness.")
    print("\n--- Step 1: Find Design Rainfall Intensity (i) ---")
    print("Iterating to balance Time of Concentration (t_c) and Rainfall Intensity (i)...")
    print(f"{'Iteration':<10} | {'Intensity i (mm/hr)':<22} | {'Concentration t_c (min)':<25}")
    print("-" * 65)

    for i in range(20):
        previous_i = design_i_mmhr
        design_tc_min = calculate_tc(design_i_mmhr, manning_n, drainage_path_L_m, cross_slope_S)
        design_i_mmhr = calculate_i(design_tc_min)
        print(f"{i+1:<10} | {design_i_mmhr:<22.2f} | {design_tc_min:<25.4f}")
        if abs(design_i_mmhr - previous_i) < tolerance:
            break

    print("\nIteration converged.")
    print(f"Final Design Rainfall Intensity (i): {design_i_mmhr:.2f} mm/hr")
    print(f"Final Time of Concentration (t_c): {design_tc_min:.4f} minutes")

    # --- Step 3: Calculate Water Film Thickness (Tw) ---
    print("\n--- Step 2: Calculate Water Film Thickness (Tw) ---")
    
    # Convert converged intensity from mm/hr to m/s for the physics formula
    i_m_per_s = design_i_mmhr / (1000 * 3600)
    
    # The formula for water film thickness in meters is:
    # Tw(m) = [ (n * i * L) / (S^0.5) ]^0.6
    # Note: q (flow) = i * L
    water_film_thickness_m = ((manning_n * i_m_per_s * drainage_path_L_m) / (cross_slope_S**0.5))**0.6
    
    # Convert final result to millimeters
    water_film_thickness_mm = water_film_thickness_m * 1000
    
    print("\nUsing the final design rainfall intensity, the water film thickness is:")
    print("Formula: Tw (mm) = 1000 * [ (n * (i_mmhr / 3600000) * L) / (S^0.5) ]^0.6\n")

    # --- Step 4: Present the Final Answer ---
    print("--- Final Equation with Numerical Values ---")
    print(
        f"Tw (mm) = 1000 * [ ({manning_n} * ({design_i_mmhr:.2f} / 3600000) * {drainage_path_L_m}) / ({cross_slope_S}**0.5) ]**0.6"
    )
    print(f"Tw (mm) = {water_film_thickness_mm:.2f}")


# Run the main function
calculate_water_film_thickness()