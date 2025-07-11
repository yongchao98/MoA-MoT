import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a roadway based on
    given parameters and standard engineering assumptions.
    """

    # Step 1: Define problem parameters and state necessary assumptions.
    # --- Given Parameters from the problem description ---
    num_lanes = 3
    lane_width_m = 3.6
    cross_slope_percent = 1.75

    # --- Assumed Parameters based on standard engineering practice ---
    # Since rainfall-duration-frequency curves are not provided, a reasonable
    # design rainfall intensity must be assumed to represent a critical storm event.
    rainfall_intensity_mmhr = 150.0  # Assumed design rainfall intensity in mm/hr

    # Manning's roughness coefficient for rough-textured asphalt pavement.
    manning_n = 0.015  # Assumed Manning's n

    print("--- Problem Parameters and Assumptions ---")
    print(f"Number of lanes per direction: {num_lanes}")
    print(f"Width of each lane: {lane_width_m} m")
    print(f"Cross-slope of the road: {cross_slope_percent}%")
    print(f"Assumed Design Rainfall Intensity (I): {rainfall_intensity_mmhr} mm/hr")
    print(f"Assumed Manning's Roughness Coefficient (n): {manning_n}\n")

    # Step 2: Calculate inputs and convert units for the formula (SI units).
    print("--- Preparing Inputs for Calculation ---")
    # Calculate the total flow path length (L), which is the width of the travelled way.
    L = num_lanes * lane_width_m
    print(f"Total Flow Path Length (L) = {num_lanes} lanes * {lane_width_m} m/lane = {L} m")

    # Convert cross-slope from percent to a decimal.
    S = cross_slope_percent / 100
    print(f"Cross-slope (S) as a decimal = {S}")

    # Convert rainfall intensity from mm/hr to m/s.
    # Conversion: 1 mm/hr = 1 / (1000 mm/m * 3600 s/hr) m/s
    I_ms = rainfall_intensity_mmhr / (1000 * 3600)
    print(f"Rainfall Intensity (I) in m/s = {I_ms:.8f} m/s\n")

    # Step 3: Apply the Kinematic Wave formula to calculate water film thickness (Tw).
    print("--- Water Film Thickness Calculation ---")
    print("The formula used is the Kinematic Wave equation for shallow sheet flow:")
    print("Tw = [ (I * L * n) / S^0.5 ] ^ 0.6\n")

    # Perform the calculation steps.
    sqrt_S = math.sqrt(S)
    numerator = I_ms * L * manning_n
    term_inside_brackets = numerator / sqrt_S
    Tw_m = math.pow(term_inside_brackets, 0.6)

    # Convert the final result from meters to millimeters.
    Tw_mm = Tw_m * 1000

    print("Substituting the values into the equation:")
    print(f"Tw (in meters) = [ ({I_ms:.8f} * {L} * {manning_n}) / ({S})**0.5 ]**0.6")
    print(f"Tw (in meters) = [ ({numerator:.8f}) / {sqrt_S:.5f} ]**0.6")
    print(f"Tw (in meters) = [ {term_inside_brackets:.8f} ]**0.6")
    print(f"Tw (in meters) = {Tw_m:.6f} m\n")

    # Step 4: Final result.
    print("--- Final Design Water Film Thickness ---")
    print(f"The calculated design water film thickness is {Tw_mm:.2f} mm.")

if __name__ == '__main__':
    calculate_water_film_thickness()
<<<4.13>>>