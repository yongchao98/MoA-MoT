import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a pavement surface.
    
    This function uses the pavement drainage equation, which is a variation of the 
    Manning's equation for shallow sheet flow. It relies on standard assumptions
    for Manning's roughness coefficient and design rainfall intensity where
    specific data is not provided.
    """
    # 1. Given and Assumed Parameters
    lane_width_m = 3.6  # Lane width in meters
    num_lanes = 3       # Number of lanes in one direction
    cross_slope_percent = 1.75  # Cross-slope in percent

    # Standard design assumption for rough-textured asphalt (FHWA HEC-22)
    manning_n = 0.016
    
    # Standard design assumption for rainfall intensity for an arterial road 
    # when hydroplaning is a concern (e.g., 10-year storm)
    rainfall_intensity_mm_hr = 150.0

    # 2. Derived Parameters and Unit Conversion
    # Total flow path length (W) in meters
    flow_path_length_W_m = num_lanes * lane_width_m
    
    # Pavement cross-slope (S) as a decimal
    cross_slope_S = cross_slope_percent / 100
    
    # Rainfall intensity (i) converted from mm/hr to m/s
    # 1 mm/hr = 1 / (1000 * 3600) m/s
    rainfall_intensity_i_m_s = rainfall_intensity_mm_hr / (1000 * 3600)

    # 3. Calculation using the Pavement Drainage Equation
    # d_w = [ (i * W * n) / S^(1/2) ]^(3/5)
    
    # Numerator of the main fraction: (i * W * n)
    numerator = rainfall_intensity_i_m_s * flow_path_length_W_m * manning_n
    
    # Denominator of the main fraction: S^(1/2)
    denominator = math.sqrt(cross_slope_S)
    
    # Water film thickness in meters (d_w)
    # The exponent 3/5 is 0.6
    d_w_m = (numerator / denominator) ** 0.6
    
    # 4. Convert result to millimeters
    d_w_mm = d_w_m * 1000

    # 5. Print the results including the final equation
    print("--- Design Parameters ---")
    print(f"Flow Path Length (W): {flow_path_length_W_m:.1f} m")
    print(f"Pavement Cross-Slope (S): {cross_slope_S}")
    print(f"Assumed Manning's n: {manning_n}")
    print(f"Assumed Design Rainfall Intensity (i): {rainfall_intensity_mm_hr} mm/hr ({rainfall_intensity_i_m_s:.2e} m/s)")
    print("\n--- Calculation ---")
    print("Formula: Water Film Thickness (mm) = [ (i * W * n) / (S^0.5) ]^0.6 * 1000")
    print("\nFinal Equation with values:")
    equation = f"Water Film Thickness (mm) = [ ({rainfall_intensity_i_m_s:.2e} m/s * {flow_path_length_W_m:.1f} m * {manning_n}) / ({cross_slope_S}^0.5) ]^0.6 * 1000"
    print(equation)
    print(f"\nDesign Water Film Thickness: {d_w_mm:.2f} mm")

if __name__ == '__main__':
    calculate_water_film_thickness()