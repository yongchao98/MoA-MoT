import math

def calculate_water_film_thickness():
    """
    Calculates the water film thickness on a pavement surface using the kinematic wave equation.
    """
    # Step 1: Define the input parameters based on the problem description.
    
    # Each lane measures 3.6 m. There are three lanes per direction.
    num_lanes = 3
    lane_width_m = 3.6
    
    # The road has a 1.75% cross-slope.
    cross_slope_percent = 1.75
    
    # Pavement is rough-textured asphalt. A typical Manning's n is assumed.
    manning_n = 0.016
    
    # A design rainfall intensity for hydroplaning analysis is assumed.
    # 100 mm/hr is a common value for critical, short-duration storms.
    rainfall_intensity_mm_hr = 100.0

    # Step 2: Calculate the terms for the kinematic wave equation in SI units.
    
    # Flow path length (L) in meters
    L_m = num_lanes * lane_width_m
    
    # Cross-slope (S) as a decimal (m/m)
    S = cross_slope_percent / 100.0
    
    # Rainfall intensity (i) in meters per second
    # Conversion: 1000 mm = 1 m, 3600 s = 1 hr
    i_m_s = rainfall_intensity_mm_hr / (1000 * 3600)

    # Step 3: Apply the kinematic wave equation: d = [ (n * L * i) / (S^0.5) ]^0.6
    
    # Calculate the numerator of the term in the brackets
    numerator = manning_n * L_m * i_m_s
    
    # Calculate the denominator of the term in the brackets
    denominator = math.sqrt(S)
    
    # Calculate the term inside the brackets
    term_in_brackets = numerator / denominator
    
    # Calculate the water depth (d) in meters
    d_m = math.pow(term_in_brackets, 0.6)
    
    # Convert the final result to millimeters
    d_mm = d_m * 1000

    # Step 4: Print the explanation and the final equation with values.
    print("Using the Kinematic Wave equation to find the water film thickness (d):")
    print("d(mm) = [ (n * L * i) / (S^0.5) ]^0.6 * 1000\n")
    print("Input Values:")
    print(f"  Manning's n (for rough asphalt): {manning_n}")
    print(f"  Flow Path Length (L): {num_lanes} lanes * {lane_width_m} m/lane = {L_m:.1f} m")
    print(f"  Cross-Slope (S): {cross_slope_percent}% = {S}")
    print(f"  Assumed Rainfall Intensity (i): {rainfall_intensity_mm_hr} mm/hr = {i_m_s:.3e} m/s\n")
    
    print("Calculation Steps:")
    final_equation_str = (
        f"d (mm) = [ ({manning_n} * {L_m:.1f} m * {i_m_s:.4e} m/s) / ({S:.4f}^0.5) ]^0.6 * 1000\n"
        f"d (mm) = [ {numerator:.3e} / {denominator:.4f} ]^0.6 * 1000\n"
        f"d (mm) = [ {term_in_brackets:.3e} ]^0.6 * 1000\n"
        f"d (mm) = {d_m:.5f} m * 1000\n"
        f"d (mm) = {d_mm:.2f} mm"
    )
    print(final_equation_str)
    
    # Return the final numeric answer for the submission format
    return round(d_mm, 2)

# Run the calculation and store the final answer
final_answer = calculate_water_film_thickness()
# The final answer is wrapped according to the required format.
# print(f"\n<<<{final_answer}>>>")