import math

def check_blackhole_entropy():
    """
    This function checks the calculation for the entropy of a supermassive black hole
    based on its angular size and distance.
    """
    # --- 1. Given parameters from the question ---
    d_pc = 10**10  # distance in parsecs
    theta_deg = 10**-17  # angular size in degrees

    # --- 2. High-precision physical constants ---
    # Using values from CODATA 2018, similar to scipy.constants
    pc_to_m = 3.085677581491367e16  # Parsec to meters
    k_B = 1.380649e-23              # Boltzmann constant (J/K)
    c = 299792458                   # Speed of light (m/s)
    G = 6.67430e-11                 # Gravitational constant (N·m²/kg²)
    hbar = 1.054571817e-34          # Reduced Planck constant (J·s)

    # --- 3. Convert parameters to SI units ---
    # Convert distance from parsecs to meters
    d_m = d_pc * pc_to_m
    
    # Convert angular size from degrees to radians
    theta_rad = theta_deg * (math.pi / 180)

    # --- 4. Calculate the radius of the event horizon (R_s) ---
    # Using the small-angle approximation, the diameter D = d * θ
    diameter = d_m * theta_rad
    
    # The radius is half the diameter
    R_s = diameter / 2

    # --- 5. Calculate the area of the event horizon (A) ---
    # The area of a sphere is A = 4 * π * R_s^2
    A = 4 * math.pi * R_s**2

    # --- 6. Calculate the Bekenstein-Hawking entropy (S) ---
    # Formula: S = (k_B * c^3 * A) / (4 * G * ħ)
    # The (4) in the numerator and denominator cancel out with A = 4*pi*R_s^2
    # S = (k_B * c^3 * (4 * math.pi * R_s**2)) / (4 * G * hbar)
    # S = (k_B * math.pi * c**3 * R_s**2) / (G * hbar)
    
    numerator = k_B * math.pi * (c**3) * (R_s**2)
    denominator = G * hbar
    entropy = numerator / denominator

    # --- 7. Check the order of magnitude ---
    # The provided answer is C) 10^62 J/K, so the expected order of magnitude is 62.
    expected_order_of_magnitude = 62
    
    # Calculate the order of magnitude of our result
    if entropy <= 0:
        return f"Incorrect. Calculated entropy is non-positive: {entropy:.4e} J/K."
        
    calculated_order_of_magnitude = math.floor(math.log10(entropy))

    # --- 8. Final verification ---
    if calculated_order_of_magnitude == expected_order_of_magnitude:
        # The options are separated by several orders of magnitude, so if the exponent
        # matches, the answer is correct.
        return "Correct"
    else:
        return (f"Incorrect. The calculated entropy is approximately {entropy:.4e} J/K. "
                f"This corresponds to an order of magnitude of 10^{calculated_order_of_magnitude}, "
                f"which does not match the answer's order of magnitude of 10^{expected_order_of_magnitude}.")

# Run the check
result = check_blackhole_entropy()
print(result)