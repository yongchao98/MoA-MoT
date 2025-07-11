import math

def calculate_angle():
    """
    Calculates the angle of particle C in the lab frame.

    The function follows these steps:
    1. Define the velocity of particle A, beta_A.
    2. Calculate the corresponding Lorentz factor, gamma_A.
    3. The relation for the tangent of the lab-frame angle theta_C is derived from Lorentz transformations:
       tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))
    4. Calculate the terms of this equation.
    5. Compute tan(theta_C) and then the angle theta_C in degrees.
    6. Print the equation with numerical values and the final result.
    """
    # Given velocity of particle A
    beta_A = 0.95

    # Calculate the Lorentz factor
    gamma_A = 1 / math.sqrt(1 - beta_A**2)

    # The term 1/sqrt(2) comes from the components of C's momentum in A's rest frame
    inv_sqrt_2 = 1 / math.sqrt(2)

    # Calculate the tangent of the angle theta_C
    tan_theta_C = inv_sqrt_2 / (gamma_A * (inv_sqrt_2 + beta_A))

    # Calculate the angle in radians
    theta_C_rad = math.atan(tan_theta_C)

    # Convert the angle to degrees
    theta_C_deg = math.degrees(theta_C_rad)

    # Print the equation with the numerical values
    print("The angle theta_C is determined by the equation:")
    print(f"tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))")
    print("\nSubstituting the numerical values:")
    print(f"tan(theta_C) = ({inv_sqrt_2:.3f}) / ({gamma_A:.3f} * ({inv_sqrt_2:.3f} + {beta_A:.3f}))")
    
    # Calculating the final value from the intermediate rounded values for display
    denominator_val = gamma_A * (inv_sqrt_2 + beta_A)
    print(f"tan(theta_C) = {inv_sqrt_2:.3f} / {denominator_val:.3f} = {tan_theta_C:.4f}")

    # Print the final result rounded to three decimal places
    print(f"\nThe angle between particle C and its parent particle A is: {theta_C_deg:.3f} degrees.")

calculate_angle()