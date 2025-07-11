import math

def solve_star_distance():
    """
    Calculates the angular distance between two stars based on their precession data.
    """
    # Given parameters
    T = 26000  # years, period of precession
    tA = -3000  # years, Star A was on the equator
    tB = 10000 # years, Star B will be on the equator
    epsilon_deg = 23.5  # degrees, Earth's axial tilt

    # Step 1: Calculate angular speed of precession
    omega_deg_per_year = 360.0 / T

    # Step 2 & 3: From the analysis of the crossing times and the coordinate swap symmetry,
    # we derive a relationship between the longitude difference (delta_lambda) and a characteristic
    # angle theta related to the star's ecliptic latitude.
    # The relation is delta_lambda = +/- 2 * theta.
    # The crossing times give another relation: delta_lambda = omega * (tA - tB) +/- 2 * theta.
    # Solving these simultaneously gives the value for theta.
    # The only consistent solution is 4 * theta = omega * (tB - tA).
    
    delta_t_eq = tB - tA
    print(f"Time difference between equator crossings: {delta_t_eq} years")

    # This leads to 4 * theta = omega * 13000 years
    theta_deg = (omega_deg_per_year * delta_t_eq) / 4.0
    print(f"The characteristic angle theta is: {theta_deg:.2f} degrees")

    # From this, delta_lambda = -2 * theta
    delta_lambda_deg = -2 * theta_deg
    print(f"The difference in ecliptic longitude (Delta Lambda) is: {delta_lambda_deg:.2f} degrees")
    
    # Step 4: Calculate the final angular distance D.
    # The formula for the distance D simplifies because cos(delta_lambda) = cos(-90) = 0.
    # cos(D) = cos^2(d), where d is the ecliptic colatitude.
    # We find d from the relation: cos(theta) = -cot(d) * cot(epsilon)

    # Convert angles to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    epsilon_rad = math.radians(epsilon_deg)

    # Calculate cot(d)
    # cot(d) = -cos(theta) / cot(epsilon) = -cos(theta) * tan(epsilon)
    cot_d = -math.cos(theta_rad) * math.tan(epsilon_rad)
    
    # Calculate cos^2(d)
    # cos^2(d) = 1 / (1 + tan^2(d)) = 1 / (1 + 1/cot^2(d)) = cot^2(d) / (1 + cot^2(d))
    cot_d_sq = cot_d**2
    cos_d_sq = cot_d_sq / (1 + cot_d_sq)

    # cos(D) = cos^2(d)
    cos_D = cos_d_sq
    print(f"The cosine of the angular distance is: {cos_D:.4f}")

    # Calculate D
    D_rad = math.acos(cos_D)
    D_deg = math.degrees(D_rad)

    # Final Answer
    print(f"\nThe equation for the final answer is arccos(tan^2({epsilon_deg}) / (tan^2({epsilon_deg}) + 2))")
    print(f"The angular distance between the two stars is {D_deg:.1f} degrees.")
    
    return D_deg

# Run the calculation
final_distance = solve_star_distance()