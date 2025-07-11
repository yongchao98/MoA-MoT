import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a roadway for hydroplaning analysis.
    """
    # Step 1: Define given variables and standard assumptions.
    
    # Given roadway characteristics
    num_lanes = 3
    lane_width_m = 3.6  # in meters
    cross_slope_percent = 1.75

    # Assumed standard design parameters
    # Manning's roughness coefficient 'n' for rough-textured asphalt pavement
    manning_n = 0.016
    # Design rainfall intensity 'i' for hydroplaning analysis (a common high-intensity value)
    rainfall_intensity_mm_hr = 100

    # Step 2: Calculate inputs for the formula in base SI units.
    
    # Drainage path length (L) in meters is the total width of the three lanes.
    L_m = num_lanes * lane_width_m

    # Cross-slope (S) as a dimensionless ratio (m/m).
    S = cross_slope_percent / 100

    # Convert rainfall intensity (i) from mm/hr to meters per second (m/s).
    # 1000 mm = 1 m, 3600 s = 1 hr
    i_m_s = rainfall_intensity_mm_hr / (1000 * 3600)

    # Step 3: Apply the kinematic wave equation for water film thickness.
    # The formula is: d = (n * L)^0.6 * i^0.6 * S^-0.3
    # The result 'd_meters' will be in meters.
    
    d_meters = (manning_n * L_m)**0.6 * (i_m_s)**0.6 * (S)**-0.3

    # Step 4: Convert the final result from meters to millimeters.
    d_mm = d_meters * 1000

    # Step 5: Print the detailed breakdown of the calculation.
    print("Calculation of Design Water Film Thickness")
    print("------------------------------------------")
    print("The calculation is based on the kinematic wave equation: d = (n * L)^0.6 * i^0.6 * S^-0.3\n")
    
    print("Input Parameters:")
    print(f"  - Number of lanes: {num_lanes}")
    print(f"  - Lane width: {lane_width_m} m")
    print(f"  - Cross-slope: {cross_slope_percent}%")
    print(f"  - Assumed Manning's n for rough asphalt: {manning_n}")
    print(f"  - Assumed design rainfall intensity: {rainfall_intensity_mm_hr} mm/hr\n")

    print("Equation with values (in base SI units):")
    print(f"  L (drainage path) = {num_lanes} * {lane_width_m} = {L_m:.2f} m")
    print(f"  S (slope) = {cross_slope_percent} / 100 = {S}")
    print(f"  i (intensity) = {rainfall_intensity_mm_hr} mm/hr = {i_m_s:.8f} m/s")
    
    print("\nFinal Equation:")
    # The user requested to output each number in the final equation.
    print(f"d (mm) = (({manning_n} * {L_m:.2f})^0.6 * ({i_m_s:.8f})^0.6 * ({S})^-0.3) * 1000")
    
    # Breaking down the terms for clarity
    term1 = (manning_n * L_m)**0.6
    term2 = i_m_s**0.6
    term3 = S**-0.3
    print(f"d (mm) = ({term1:.4f} * {term2:.6f} * {term3:.4f}) * 1000")
    print(f"d (mm) = {d_meters:.6f} * 1000")

    print("\nResult:")
    print(f"The design water film thickness is {d_mm:.2f} mm.")

# Execute the function
calculate_water_film_thickness()