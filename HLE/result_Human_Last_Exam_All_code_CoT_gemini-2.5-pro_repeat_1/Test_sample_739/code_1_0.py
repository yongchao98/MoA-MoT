import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a pavement surface
    based on road geometry and assumed rainfall conditions.
    """

    # Step 1: Define the given and assumed parameters.
    # Number of lanes in one direction
    N = 3
    # Width of each lane in meters
    W_L = 3.6  # m
    # Cross-slope of the pavement in percent
    S_percent = 1.75  # %
    
    # Manning's roughness coefficient for rough-textured asphalt pavement (assumed)
    n = 0.016  # s/m^(1/3)
    
    # Design rainfall intensity (assumed for a 10-year, 5-minute storm)
    i_mmhr = 200  # mm/hr

    # Step 2: Calculate intermediate values and convert units for the formula.
    # Total drainage path length (L) in meters
    L = N * W_L
    
    # Convert cross-slope from percent to a decimal value (m/m)
    S = S_percent / 100
    
    # Convert rainfall intensity from mm/hr to m/s for use in the SI unit-based formula
    # 1 hr = 3600 s, 1 m = 1000 mm
    i_ms = i_mmhr / (1000 * 3600)
    
    # Step 3: Apply the Kinematic Wave equation for water film thickness.
    # The formula is d = [ (n * i * L) / sqrt(S) ]^(3/5), where d is in meters.
    
    # Calculate the term inside the brackets
    numerator = n * i_ms * L
    denominator = math.sqrt(S)
    term_inside = numerator / denominator
    
    # The exponent is 3/5
    exponent = 3.0 / 5.0
    
    # Calculate the water film thickness in meters
    d_m = math.pow(term_inside, exponent)
    
    # Convert the water film thickness to millimeters
    d_mm = d_m * 1000
    
    # Step 4: Print the inputs, the equation, and the final result.
    print("To determine the design water film thickness, we use the Kinematic Wave equation:")
    print("d_mm = [ (n * i_ms * L) / S^(1/2) ]^(3/5) * 1000")
    print("\n--- Input Parameters ---")
    print(f"Number of lanes per direction: {N}")
    print(f"Lane width: {W_L} m")
    print(f"Cross-slope: {S_percent}%")
    print(f"Assumed Manning's n: {n}")
    print(f"Assumed rainfall intensity: {i_mmhr} mm/hr")

    print("\n--- Calculation using the formula ---")
    print(f"1. Drainage Path Length (L) = {N} * {W_L} = {L:.1f} m")
    print(f"2. Cross-Slope (S) = {S_percent} / 100 = {S:.4f}")
    print(f"3. Rainfall Intensity (i) = {i_mmhr} mm/hr = {i_ms:.7f} m/s")
    
    print("\nFinal Equation with values:")
    # Using scientific notation for clarity in the equation string
    print(f"d_mm = [ ({n} * {i_ms:.2e} * {L:.1f}) / ({S:.4f}^0.5) ]^{exponent:.1f} * 1000")
    print(f"d_mm = [ ({numerator:.2e}) / {denominator:.4f} ]^{exponent:.1f} * 1000")
    print(f"d_mm = [ {term_inside:.2e} ]^{exponent:.1f} * 1000")
    print(f"d_mm = {d_m:.5f} m * 1000")
    
    print(f"\nResult:")
    print(f"The design water film thickness is {d_mm:.2f} mm.")

    return d_mm

# Execute the function and store the final answer.
final_answer = calculate_water_film_thickness()
# The final answer is wrapped in <<<>>> as requested.
# The print statements within the function will be shown to the user.
# This last line with the final answer will not be printed to the user console.
final_answer_formatted = f"<<<{final_answer:.2f}>>>"

# The prompt asks for a single code block, so this comment and the final_answer_formatted line
# are for thought process. The code block will end before this.

if __name__ == '__main__':
    calculate_water_film_thickness()