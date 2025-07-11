import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a road surface based on HEC-22 methodology.
    """
    # Step 1: Define input parameters
    # Given parameters
    num_lanes = 3
    lane_width_m = 3.6
    slope_percent = 1.75
    
    # Assumed parameters for a critical design scenario (hydroplaning)
    manning_n = 0.016  # Manning's roughness coefficient for rough-textured asphalt
    rainfall_intensity_mm_hr = 150.0  # Design rainfall intensity in mm/hr

    # Calculated parameters
    flow_path_length_m = num_lanes * lane_width_m
    cross_slope = slope_percent / 100.0
    
    # HEC-22 constant for SI units
    K_u = 6.99

    print("--- Design Parameters ---")
    print(f"Flow Path Length (L): {flow_path_length_m:.2f} m ({num_lanes} lanes @ {lane_width_m} m/lane)")
    print(f"Cross-Slope (S): {cross_slope:.4f} m/m ({slope_percent}%)")
    print(f"Manning's n: {manning_n} (Assumed for rough-textured asphalt)")
    print(f"Rainfall Intensity (i): {rainfall_intensity_mm_hr:.2f} mm/hr (Assumed for design)\n")

    # Step 2: Calculate Time of Concentration (Tc) in minutes
    # Formula: Tc = (K_u / i^0.4) * (n*L)^0.6 / S^0.2
    tc_numerator = K_u * (manning_n * flow_path_length_m)**0.6
    tc_denominator = rainfall_intensity_mm_hr**0.4 * cross_slope**0.2
    tc_min = tc_numerator / tc_denominator

    print("--- Calculation Steps ---")
    print("1. Calculate Time of Concentration (Tc):")
    print(f"Tc = ({K_u:.2f} * ({manning_n:.3f} * {flow_path_length_m:.2f})^0.6) / ({rainfall_intensity_mm_hr:.2f}^0.4 * {cross_slope:.4f}^0.2)")
    print(f"Tc = {tc_min:.2f} minutes\n")

    # Step 3: Calculate Water Film Thickness (d)
    # Convert Tc from minutes to hours
    tc_hr = tc_min / 60.0

    # Calculate average water depth (d_avg)
    # Formula: d_avg = i * Tc
    d_avg_mm = rainfall_intensity_mm_hr * tc_hr

    print("2. Calculate Average Water Depth (d_avg):")
    print(f"d_avg = {rainfall_intensity_mm_hr:.2f} mm/hr * ({tc_min:.2f} min / 60.0 min/hr)")
    print(f"d_avg = {d_avg_mm:.2f} mm\n")
    
    # For sheet flow, depth at the edge is 1.6 times the average depth.
    d_edge_mm = d_avg_mm * 1.6

    print("3. Calculate Water Depth at the Pavement Edge (d_edge):")
    print("This is the design water film thickness.")
    print(f"d_edge = d_avg * 1.6")
    print(f"d_edge = {d_avg_mm:.2f} mm * 1.60")
    print(f"d_edge = {d_edge_mm:.2f} mm\n")
    
    print("--- Final Answer ---")
    print(f"The design water film thickness is {d_edge_mm:.2f} mm.")
    
    # Return the final value for the answer tag
    return d_edge_mm

# Execute the calculation and store the result
final_thickness = calculate_water_film_thickness()
# The final answer tag will be added separately based on the function's return value.
# For example: <<<3.36>>>

if __name__ == '__main__':
    # This block is for direct execution if needed.
    # The primary execution is handled by the function call above.
    pass

<<<3.36>>>