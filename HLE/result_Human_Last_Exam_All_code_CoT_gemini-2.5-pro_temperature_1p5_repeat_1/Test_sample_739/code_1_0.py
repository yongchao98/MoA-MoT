import math

def calculate_water_film_thickness():
    """
    Calculates the water film thickness on a road surface based on road design
    parameters and assumed rainfall, using the Kinematic Wave equation.
    """
    # --- Step 1: Define Given and Assumed Parameters ---

    # Manning's roughness coefficient (n) for rough-textured asphalt pavement
    n = 0.016

    # Number of lanes in one direction
    num_lanes = 3

    # Lane width in meters
    lane_width = 3.6

    # Cross-slope in percent
    cross_slope_percent = 1.75

    # Design rainfall intensity (i) in mm/hr. This is a standard assumption
    # for hydroplaning analysis as it was not provided in the prompt.
    rainfall_intensity_mm_hr = 100

    # --- Step 2: Prepare Values for the Formula ---

    # Calculate the total flow path length (L) in meters
    L = num_lanes * lane_width

    # Convert cross-slope from percent to a decimal value (S)
    S = cross_slope_percent / 100

    # Convert rainfall intensity from mm/hr to m/s (i_ms)
    # 1 mm/hr = 1 / (1000 * 3600) m/s
    i_ms = rainfall_intensity_mm_hr / (1000 * 3600)

    # --- Step 3: Calculate the Water Film Thickness ---

    # The formula for water film thickness (Td) in mm is:
    # Td = 1000 * [ (n * L * i_ms) / (S^0.5) ]^0.6
    
    # Calculate intermediate terms for clarity
    numerator = n * L * i_ms
    denominator = math.sqrt(S)
    
    # Calculate depth in meters first
    depth_m = (numerator / denominator) ** 0.6
    
    # Convert the final result to millimeters (mm)
    Td_mm = depth_m * 1000

    # --- Step 4: Output the Results ---
    print("--- Design Parameters ---")
    print(f"Manning's Roughness Coefficient (n): {n}")
    print(f"Total Flow Path Length (L): {L:.1f} m ({num_lanes} lanes * {lane_width} m/lane)")
    print(f"Cross-Slope (S): {S} (from {cross_slope_percent}%)")
    print(f"Assumed Rainfall Intensity (i): {rainfall_intensity_mm_hr} mm/hr")
    print("\n--- Calculation Steps ---")
    print("The final equation for water film thickness (Td) with the given values is:")
    
    # Print the equation with all numbers substituted
    print(f"\nTd(mm) = 1000 * [({n} * {L:.1f} * {i_ms:.8f}) / ({S}**0.5)]**0.6")
    
    # Print the equation with the numerator and denominator calculated
    print(f"Td(mm) = 1000 * [({numerator:.8f}) / ({denominator:.5f}) ]**0.6")

    # Print the equation with the fraction calculated
    print(f"Td(mm) = 1000 * [{(numerator/denominator):.8f}]**0.6")
    
    # Print the equation after applying the exponent
    print(f"Td(mm) = 1000 * {depth_m:.6f}")
    
    print("\n--- Final Result ---")
    print(f"The design water film thickness is: {Td_mm:.2f} mm")


if __name__ == "__main__":
    calculate_water_film_thickness()