import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a road surface based on given parameters.
    """
    # 1. Define input parameters based on the problem description.
    num_lanes = 3
    lane_width_m = 3.6
    cross_slope_percent = 1.75

    # 2. Determine derived parameters and make necessary engineering assumptions.
    # Flow path length (L) in meters
    L = num_lanes * lane_width_m
    
    # Cross-slope (S) in decimal form (m/m)
    S = cross_slope_percent / 100
    
    # Manning's roughness coefficient (n) for rough-textured asphalt.
    n = 0.016
    
    # Rainfall intensity (i) in mm/hr.
    # This is a critical assumption as the IDF curve was not provided.
    # We assume a 10-year, 5-minute storm with an intensity of 150 mm/hr.
    i_mmhr = 150
    
    # 3. Convert units for the formula.
    # Convert rainfall intensity from mm/hr to m/s.
    i_ms = i_mmhr / (1000 * 3600)
    
    # 4. Calculate the water film thickness using the kinematic wave formula.
    # The formula is: d = [ (i * L * n) / sqrt(S) ] ^ 0.6
    try:
        term_in_brackets = (i_ms * L * n) / math.sqrt(S)
        d_m = term_in_brackets ** 0.6
        d_mm = d_m * 1000
    except ValueError as e:
        print(f"Error during calculation: {e}")
        return

    # 5. Print the results, including the final equation with values.
    print("This script calculates the design water film thickness for a roadway.")
    print("------------------------------------------------------------------")
    print("Given and Assumed Parameters:")
    print(f"  Flow Path Length (L): {L:.1f} m (3 lanes at 3.6 m each)")
    print(f"  Cross-Slope (S): {S:.4f} m/m ({cross_slope_percent}%)")
    print(f"  Manning's Roughness (n): {n:.3f} (for rough asphalt)")
    print(f"  Assumed Rainfall Intensity (i): {i_mmhr} mm/hr")
    print("------------------------------------------------------------------")
    
    print("\nThe water film thickness (d) is calculated using the formula:")
    print("d_meters = [ (i_ms * L * n) / sqrt(S) ] ^ 0.6\n")
    
    print("Substituting the values into the equation:")
    # We use f-strings to format the numbers neatly within the equation string.
    equation = (
        f"d_mm = [ ({i_ms:.8f} * {L:.1f} * {n:.3f}) / sqrt({S:.4f}) ]^0.6 * 1000"
    )
    print(equation)

    inner_calc1 = i_ms * L * n
    inner_calc2 = math.sqrt(S)
    final_calc_step = (inner_calc1 / inner_calc2)**0.6 * 1000
    
    print(f"d_mm = [ ({inner_calc1:.8f}) / {inner_calc2:.5f} ]^0.6 * 1000")
    print(f"d_mm = [ {(inner_calc1 / inner_calc2):.8f} ]^0.6 * 1000")
    print(f"d_mm = {final_calc_step:.2f} mm\n")
    
    print("Final Answer:")
    print(f"The design water film thickness is {d_mm:.2f} mm.")
    print(f"<<<{d_mm:.2f}>>>")

if __name__ == "__main__":
    calculate_water_film_thickness()
