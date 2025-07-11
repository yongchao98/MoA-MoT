import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a road surface based on
    given and assumed parameters.
    """
    # Step 1: Define problem parameters and assumptions
    # Given parameters
    num_lanes = 3
    lane_width_m = 3.6
    cross_slope_percent = 1.75

    # Assumed parameters based on standard engineering practice
    # Manning's roughness coefficient 'n' for rough-textured asphalt pavement
    n = 0.016
    # Design rainfall intensity 'i' in mm/hr for hydroplaning analysis
    i_mm_hr = 100.0

    # Constant 'K' for the formula in metric units
    K = 0.0003

    # Step 2: Prepare values for the formula
    # Calculate the total flow path length 'L' in meters
    L_m = num_lanes * lane_width_m
    # Convert the cross-slope 'S' from percent to a decimal
    S = cross_slope_percent / 100.0

    # Step 3: Apply the water film thickness formula
    # Formula: d = K * (n * L * i)^(3/5) / S^(3/10)
    # where:
    # d = water film thickness (mm)
    # K = metric constant (0.0003)
    # n = Manning's roughness coefficient
    # L = Flow path length (m)
    # i = Rainfall intensity (mm/hr)
    # S = Cross-slope (m/m)

    numerator_base = n * L_m * i_mm_hr
    denominator_base = S
    
    water_film_thickness_mm = (K * (numerator_base ** (3/5))) / (denominator_base ** (3/10))

    # Step 4: Print the detailed calculation and result
    print("Calculation of Design Water Film Thickness")
    print("-" * 40)
    print(f"Inputs and Assumptions:")
    print(f"  - Flow Path Length (L): {num_lanes} lanes * {lane_width_m} m/lane = {L_m:.2f} m")
    print(f"  - Cross-Slope (S): {cross_slope_percent}% = {S}")
    print(f"  - Manning's n (Rough Asphalt): {n}")
    print(f"  - Rainfall Intensity (i): {i_mm_hr} mm/hr (Assumed)")
    print(f"  - Constant (K): {K}\n")
    
    print("Applying the formula: d = K * (n * L * i)^(3/5) / S^(3/10)\n")
    
    # Print the equation with all the numerical values substituted
    print("Final Equation:")
    print(f"d = {K} * ({n} * {L_m:.2f} * {i_mm_hr})^(3/5) / {S}^(3/10)")
    
    # Calculate intermediate values for clarity in the final printout
    term1 = numerator_base ** (3/5)
    term2 = denominator_base ** (3/10)
    print(f"d = {K} * {term1:.4f} / {term2:.4f}")
    
    final_numerator = K * term1
    print(f"d = {final_numerator:.6f} / {term2:.4f}")
    
    print(f"\nResult:")
    print(f"The design water film thickness is {water_film_thickness_mm:.2f} mm.")

if __name__ == "__main__":
    calculate_water_film_thickness()
<<<5.15>>>