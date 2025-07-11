import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a roadway based on given parameters.
    """
    # Step 1: Define the given parameters and state assumptions.
    lane_width_m = 3.6
    num_lanes = 3
    cross_slope_percent = 1.75
    
    # Manning's roughness coefficient (n) for rough-textured asphalt pavement.
    # This is a standard value used in hydraulic calculations for such surfaces.
    manning_n = 0.015
    
    # Assumed design rainfall intensity (I) in mm/hr.
    # This value is typically determined from Intensity-Duration-Frequency (IDF) curves.
    # Since they are not provided, a common value for hydroplaning analysis is assumed.
    rainfall_intensity_mm_hr = 150.0

    print("--- Design Water Film Thickness Calculation ---")
    print("\nStep 1: Define Input Parameters and Assumptions")
    print(f"Number of lanes per direction = {num_lanes}")
    print(f"Lane width = {lane_width_m} m")
    print(f"Cross-slope = {cross_slope_percent}%")
    print(f"Pavement type: Rough-textured asphalt")
    print(f"Assumed Manning's roughness coefficient (n) = {manning_n}")
    print(f"Assumed design rainfall intensity (I) = {rainfall_intensity_mm_hr} mm/hr")

    # Step 2: Calculate the parameters for the water film thickness formula.
    # Calculate Flow Path Length (L) in meters.
    L = num_lanes * lane_width_m
    # Calculate Cross-slope (S) in m/m.
    S = cross_slope_percent / 100.0
    # Convert Rainfall Intensity (i) from mm/hr to m/s.
    i_m_per_s = rainfall_intensity_mm_hr / (1000 * 3600)

    print("\nStep 2: Prepare Parameters for the Calculation Formula")
    print("The formula for water film thickness (d_w) is: d_w = [ (n * L * i) / sqrt(S) ]^(3/5)")
    print("\nCalculated components for the formula:")
    print(f"Flow Path Length (L) = {num_lanes} lanes * {lane_width_m} m/lane = {L:.2f} m")
    print(f"Cross-slope (S) = {cross_slope_percent}% / 100 = {S}")
    print(f"Rainfall Intensity (i) = {rainfall_intensity_mm_hr} mm/hr / (1000 mm/m * 3600 s/hr) = {i_m_per_s:.8f} m/s")

    # Step 3: Calculate the water film thickness.
    term_in_parenthesis_numerator = manning_n * L * i_m_per_s
    term_in_parenthesis_denominator = S**0.5
    term_in_parenthesis = term_in_parenthesis_numerator / term_in_parenthesis_denominator
    d_w_m = term_in_parenthesis**(3/5)
    d_w_mm = d_w_m * 1000

    print("\nStep 3: Substitute Values into the Equation")
    print(f"d_w (m) = [ ({manning_n} * {L:.2f} * {i_m_per_s:.8f}) / ({S}**0.5) ]^(3/5)")
    print(f"d_w (m) = [ {term_in_parenthesis_numerator:.8f} / {term_in_parenthesis_denominator:.6f} ]^(3/5)")
    print(f"d_w (m) = [ {term_in_parenthesis:.8f} ]^(3/5)")
    print(f"d_w (m) = {d_w_m:.6f} m")

    # Step 4: Final result in mm.
    print("\nStep 4: Determine Final Answer in Millimeters")
    print(f"Design Water Film Thickness = {d_w_m:.6f} m * 1000 mm/m")
    print(f"Design Water Film Thickness = {d_w_mm:.2f} mm")

if __name__ == '__main__':
    calculate_water_film_thickness()