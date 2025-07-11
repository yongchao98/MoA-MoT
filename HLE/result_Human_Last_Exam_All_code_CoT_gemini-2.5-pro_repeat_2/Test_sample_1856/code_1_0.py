import math

def calculate_lab_frame_angle():
    """
    Calculates the angle of particle C in the lab frame based on its emission
    in the rest frame of a moving parent particle A.
    """
    # Given values
    beta_A = 0.95
    # The emission angle in the rest frame of A (theta*) is 45 degrees,
    # as derived from the given four-momentum C_mu*.
    theta_star_deg = 45.0
    theta_star_rad = math.radians(theta_star_deg)

    # Step 1: Calculate the Lorentz factor, gamma
    try:
        gamma_A = 1 / math.sqrt(1 - beta_A**2)
    except ValueError:
        print("Error: beta_A must be less than 1.")
        return

    # Step 2: Use the angle transformation formula for an ultra-relativistic particle
    # tan(theta_C) = sin(theta_star) / (gamma * (cos(theta_star) + beta))
    sin_theta_star = math.sin(theta_star_rad)
    cos_theta_star = math.cos(theta_star_rad)

    # Calculate tan(theta_C) by breaking down the equation
    term_in_parentheses = cos_theta_star + beta_A
    denominator = gamma_A * term_in_parentheses
    tan_theta_C = sin_theta_star / denominator

    # Step 3: Calculate the final angle in degrees
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)

    # Print the detailed calculation process
    print("This script calculates the emission angle of particle C in the lab frame.")
    print("\nGiven values:")
    print(f"Velocity of parent particle A (β_A): {beta_A}")
    print(f"Emission angle in A's rest frame (θ*): {theta_star_deg}°")
    
    print("\nStep 1: Calculate the Lorentz factor (γ_A)")
    print(f"γ_A = 1 / sqrt(1 - {beta_A}**2)")
    print(f"γ_A = {gamma_A:.4f}")

    print("\nStep 2: Calculate tan(θ_C) using the formula tan(θ_C) = sin(θ*) / (γ_A * (cos(θ*) + β_A))")
    print(f"Equation: tan(θ_C) = sin({theta_star_deg}°) / ({gamma_A:.4f} * (cos({theta_star_deg}°) + {beta_A}))")
    print(f"Substituting values: tan(θ_C) = {sin_theta_star:.4f} / ({gamma_A:.4f} * ({cos_theta_star:.4f} + {beta_A}))")
    print(f"tan(θ_C) = {sin_theta_star:.4f} / ({gamma_A:.4f} * {term_in_parentheses:.4f})")
    print(f"tan(θ_C) = {sin_theta_star:.4f} / {denominator:.4f}")
    print(f"tan(θ_C) = {tan_theta_C:.6f}")

    print("\nStep 3: Calculate the angle θ_C in degrees")
    print(f"θ_C = arctan({tan_theta_C:.6f})")
    print(f"\nThe final angle of particle C in the lab frame is: {theta_C_deg:.3f}°")

if __name__ == "__main__":
    calculate_lab_frame_angle()