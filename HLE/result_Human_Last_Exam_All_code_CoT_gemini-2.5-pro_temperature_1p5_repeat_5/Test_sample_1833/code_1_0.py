import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given parameters
    m = 1.0  # kg
    M = 1.0  # kg
    g = 9.8  # m/s^2
    theta_deg = 60.0

    # Convert angle from degrees to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Calculate trigonometric values needed for the formula
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    sin_sq_theta = sin_theta**2
    cos_sq_theta = cos_theta**2

    # Step 1: Calculate the L*omega^2 term (let's call it L_omega_sq)
    # This comes from the conservation of energy principle.
    # For the general case: Lω² = (2*g*sinθ) / [ (m/(m+M))sin²θ + cos²θ ]
    # For our case m=M=1: m/(m+M) = 0.5
    l_omega_sq_numerator = 2 * g * sin_theta
    l_omega_sq_denominator = (m / (m + M)) * sin_sq_theta + cos_sq_theta
    l_omega_sq = l_omega_sq_numerator / l_omega_sq_denominator

    # Step 2: Calculate the tension T using the derived formula for the m=M case
    # T = (M*g*sinθ + M*Lω²) / (1 + cos²θ)
    tension_numerator = M * g * sin_theta + M * l_omega_sq
    tension_denominator = 1 + cos_sq_theta
    tension = tension_numerator / tension_denominator

    # Output the components of the final calculation
    print("--- Calculation Breakdown ---")
    print(f"Angle (θ): {theta_deg}°")
    print(f"M (mass of object): {M:.2f} kg")
    print(f"g (acceleration due to gravity): {g:.2f} m/s^2")
    print(f"sin(θ): {sin_theta:.4f}")
    print(f"cos²(θ): {cos_sq_theta:.4f}")
    print(f"Lω²: {l_omega_sq:.4f} m/s²")
    print("\nFinal Equation: T = (M * g * sin(θ) + M * Lω²) / (1 + cos²(θ))")
    print(f"T = ({M:.2f} * {g:.2f} * {sin_theta:.4f} + {M:.2f} * {l_omega_sq:.4f}) / (1 + {cos_sq_theta:.4f})")
    print(f"T = ({tension_numerator:.4f}) / ({tension_denominator:.4f})")
    
    # Print the final answer rounded to two decimal places
    print("\n--- Final Answer ---")
    print(f"The tension in the string is: {tension:.2f} Newtons")

calculate_tension()
<<<28.52>>>