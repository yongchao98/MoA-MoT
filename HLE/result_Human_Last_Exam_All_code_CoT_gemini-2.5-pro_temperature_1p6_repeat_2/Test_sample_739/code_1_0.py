import math

# --- Parameters ---

# Road geometry from the problem description
lane_width_m = 3.6
num_lanes = 3
cross_slope_percent = 1.75

# Assumed parameters based on standard engineering practice
# Manning's roughness coefficient 'n' for rough-textured asphalt
n = 0.016
# Design rainfall intensity 'i' in mm/hr. A high value is assumed for hydroplaning analysis
# as the Rainfall-Duration-Frequency curves were not provided.
i_mm_hr = 100.0


# --- Calculations ---

# 1. Calculate the total flow path length (L) across the lanes
L_m = lane_width_m * num_lanes

# 2. Convert the cross-slope (S) from percent to a decimal
S_decimal = cross_slope_percent / 100

# 3. Convert the rainfall intensity (i) from mm/hr to m/s for the formula
# 1 hour = 3600 seconds, 1 meter = 1000 millimeters
i_m_s = i_mm_hr / (1000 * 3600)

# 4. Apply the water film thickness formula: d = [ (n * L * i) / (S^0.5) ]^(3/5)
# where d is in meters
try:
    numerator = n * L_m * i_m_s
    denominator = math.sqrt(S_decimal)
    base = numerator / denominator
    exponent = 3.0 / 5.0
    d_meters = math.pow(base, exponent)

    # 5. Convert the final depth from meters to millimeters
    d_mm = d_meters * 1000

    # --- Output ---
    print("Determining Design Water Film Thickness\n")
    print(f"--- Input Parameters ---")
    print(f"Flow Path Length (L): {L_m:.2f} m")
    print(f"Cross-Slope (S): {S_decimal:.4f} m/m")
    print(f"Manning's n (assumed): {n}")
    print(f"Rainfall Intensity (i, assumed): {i_mm_hr} mm/hr\n")

    print(f"--- Calculation ---")
    print("Using the formula: d_mm = [ (n * L * i_m_s) / (S^0.5) ]^(3/5) * 1000")
    
    # Printing the equation with all the numerical values substituted
    equation_str = (
        f"d_mm = [ ({n} * {L_m:.2f} * {i_m_s:.5f}) / ({S_decimal:.4f}^0.5) ]^({exponent}) * 1000"
    )
    print("\nFinal Equation:")
    print(equation_str)

    print(f"\n--- Result ---")
    print(f"The calculated design water film thickness is: {d_mm:.2f} mm")

except (ValueError, ZeroDivisionError) as e:
    print(f"An error occurred during calculation: {e}")
    print("Please check the input parameters, especially the cross-slope which cannot be zero.")
