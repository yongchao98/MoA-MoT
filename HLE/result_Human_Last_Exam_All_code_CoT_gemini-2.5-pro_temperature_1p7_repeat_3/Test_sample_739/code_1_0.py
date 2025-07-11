import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a road surface based on
    its geometry, material, and a design rainfall intensity.
    """
    # Step 1: Define Given and Assumed Parameters
    # --- Road Geometry ---
    # Number of lanes in one direction
    num_lanes = 3
    # Width of each lane in meters
    lane_width_m = 3.6
    # Cross-slope of the pavement as a percentage
    cross_slope_percent = 1.75

    # --- Pavement & Rainfall Properties ---
    # Manning's roughness coefficient 'n' for rough-textured asphalt (a standard assumed value)
    manning_n = 0.016
    # Assumed design rainfall intensity in mm/hr, as rainfall curves were not provided.
    # 150 mm/hr is a common value for a 10-year storm in many urban areas.
    i_mm_per_hr = 150.0

    # Step 2: Derive Variables for the Calculation
    # Calculate the total flow path length (L) in meters
    L_m = num_lanes * lane_width_m

    # Convert cross-slope (S) from a percentage to a decimal value
    S = cross_slope_percent / 100.0

    # Convert rainfall intensity (i) from mm/hr to m/s for use in the SI-unit based formula
    # 1 hour = 3600 seconds, 1 meter = 1000 millimeters
    i_m_per_s = i_mm_per_hr / (1000 * 3600)

    # Step 3: Apply the Manning's Equation for Shallow Flow
    # The formula for water film thickness (d) is: d = ((i * L * n) / S^0.5)^(3/5)
    # where 3/5 is 0.6 and S^0.5 is the square root of S.

    # Calculate the numerator of the inner term: i * L * n
    numerator = i_m_per_s * L_m * manning_n

    # Calculate the denominator of the inner term: sqrt(S)
    denominator = math.sqrt(S)

    # Calculate the water film thickness in meters
    d_m = (numerator / denominator) ** (3.0 / 5.0)

    # Convert the final water film thickness from meters to millimeters for the result
    d_mm = d_m * 1000

    # Step 4: Print the detailed calculation and the final result
    print("--- Water Film Thickness Calculation ---")
    print("\n[1] Input Parameters:")
    print(f"  - Flow Path Length (L): {num_lanes} lanes * {lane_width_m} m/lane = {L_m:.1f} m")
    print(f"  - Pavement Cross-Slope (S): {cross_slope_percent}% = {S}")
    print(f"  - Manning's n (Rough Asphalt): {manning_n}")
    print(f"  - Assumed Rainfall Intensity (i): {i_mm_per_hr} mm/hr = {i_m_per_s:.2e} m/s")

    print("\n[2] Formula Used:")
    print("  d = ((i * L * n) / (S^0.5)) ^ (3/5)")

    print("\n[3] Calculation with Values:")
    print(f"  d_m = (({i_m_per_s:.4e} * {L_m:.1f} * {manning_n}) / ({S}**0.5))**(3/5)")
    print(f"  d_m = (({numerator:.4e}) / ({denominator:.4f}))**0.6")
    print(f"  d_m = {d_m:.5f} m")

    print("\n[4] Final Answer:")
    print("The design water film thickness at the outer edge of the pavement is:")
    print(f"{d_mm:.2f} mm")


# Execute the function
if __name__ == "__main__":
    calculate_water_film_thickness()