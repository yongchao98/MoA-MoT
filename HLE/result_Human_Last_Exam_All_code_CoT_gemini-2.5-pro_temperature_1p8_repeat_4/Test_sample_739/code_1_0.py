import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a pavement surface.
    """

    # 1. Define parameters based on the problem statement and engineering standards.

    # Drainage Path Length (L) in meters
    # The road has 3 lanes sloping outward, each 3.6 m wide.
    num_lanes = 3
    lane_width = 3.6  # in meters
    L = num_lanes * lane_width

    # Pavement Cross-Slope (S) as a decimal value (m/m)
    cross_slope_percent = 1.75
    S = cross_slope_percent / 100

    # Manning's Roughness Coefficient (n)
    # A standard value for "rough-textured asphalt pavement" is used.
    n = 0.015

    # Rainfall Intensity (i) in mm/hr
    # Since this is not provided, a common design value for hydroplaning
    # analysis (short-duration, high-intensity storm) is assumed.
    i = 150  # in mm/hr

    # Formula Coefficient (C) for metric units
    # This coefficient is for the Anderson/Gallaway empirical formula when using metric units.
    C = 0.0135

    # 2. Calculate the water film thickness (Td) using the formula:
    # Td = C * (n * L)^0.6 * i^0.4 * S^-0.2
    td = C * (n * L)**0.6 * i**0.4 * S**-0.2

    # 3. Print the calculation steps and the final result.
    print("Calculating the design water film thickness (Td) using an empirical formula for hydroplaning analysis.\n")
    print("Formula: Td = C * (n * L)^0.6 * i^0.4 * S^-0.2\n")

    print("Defined Parameters:")
    print(f"- Coefficient for Metric Units (C): {C}")
    print(f"- Manning's n for Rough-Textured Asphalt (n): {n}")
    print(f"- Drainage Path Length (L): {L} m")
    print(f"- Pavement Cross-Slope (S): {S} m/m")
    print(f"- Assumed Design Rainfall Intensity (i): {i} mm/hr\n")

    print("Final Equation with substituted values:")
    print(f"Td (mm) = {C} * ({n} * {L})^0.6 * {i}^0.4 * {S}^-0.2\n")
    
    # Calculate intermediate terms to show the breakdown
    term1_val = n * L
    term1_res = term1_val**0.6
    term2_res = i**0.4
    term3_res = S**-0.2
    
    print("Breaking down the calculation:")
    print(f"({n} * {L})^0.6 = ({term1_val:.3f})^0.6 = {term1_res:.4f}")
    print(f"{i}^0.4 = {term2_res:.4f}")
    print(f"{S}^-0.2 = {term3_res:.4f}\n")

    print(f"Td (mm) = {C} * {term1_res:.4f} * {term2_res:.4f} * {term3_res:.4f}")
    print(f"The final calculated design water film thickness is: {td:.2f} mm")


# Execute the calculation
calculate_water_film_thickness()
<<<8.47>>>