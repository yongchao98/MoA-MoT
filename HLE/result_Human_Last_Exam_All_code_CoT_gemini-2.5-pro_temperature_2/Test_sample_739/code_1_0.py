import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a road surface using the Kinematic Wave equation.
    """
    # Step 1 & 2: Define given parameters and calculate drainage path length
    lanes_per_direction = 3
    lane_width_m = 3.6
    cross_slope_percent = 1.75
    
    # The drainage path is the total width of the lanes from the median to the outer edge.
    L_m = lanes_per_direction * lane_width_m
    
    # The cross-slope as a decimal.
    S = cross_slope_percent / 100

    # Step 3: Define assumed parameters based on standard engineering practice.
    # Manning's n for rough-textured asphalt.
    n = 0.015
    # Assumed design rainfall intensity in mm/hr (typical for a 10-year storm).
    rainfall_intensity_mm_hr = 150.0

    # Convert rainfall intensity from mm/hr to m/s for the formula.
    # 1 hr = 3600 s, 1 m = 1000 mm
    i_m_s = rainfall_intensity_mm_hr / (1000 * 3600)

    # Step 4: Apply the Kinematic Wave formula.
    # d = [ (n * i * L) / (S^0.5) ] ^ 0.6
    
    term_in_brackets = (n * i_m_s * L_m) / math.sqrt(S)
    d_m = term_in_brackets ** 0.6
    
    # Step 5: Convert the result from meters to millimeters.
    d_mm = d_m * 1000

    # Print the inputs and the final result clearly.
    print("--- Design Parameters ---")
    print(f"Number of lanes (one direction): {lanes_per_direction}")
    print(f"Lane width: {lane_width_m} m")
    print(f"Total drainage path length (L): {L_m:.1f} m")
    print(f"Pavement cross-slope (S): {cross_slope_percent}% or {S}")
    print(f"Manning's roughness coefficient (n): {n} (assumed for rough asphalt)")
    print(f"Design rainfall intensity (i): {rainfall_intensity_mm_hr} mm/hr or {i_m_s:.6f} m/s (assumed)")
    print("\n--- Calculation ---")
    print("The formula for water film thickness (d) is: d = [ (n * i * L) / (S^0.5) ] ^ 0.6")
    print("\nPlugging in the values:")
    print(f"d = [ ({n} * {i_m_s:.6f} * {L_m:.1f}) / ({S}^0.5) ] ^ 0.6")
    print(f"d = [ ({n * i_m_s * L_m:.8f}) / ({math.sqrt(S):.5f}) ] ^ 0.6")
    print(f"d = [ {term_in_brackets:.8f} ] ^ 0.6")
    print(f"d = {d_m:.6f} m")
    
    print("\n--- Final Answer ---")
    print(f"The design water film thickness is: {d_mm:.2f} mm")

# Run the calculation
calculate_water_film_thickness()
<<<4.45>>>