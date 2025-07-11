import math

def calculate_angle():
    """
    Calculates the angle between particle C and its parent particle A in the lab frame.
    """
    # Given velocity of particle A in natural units.
    beta_A = 0.95

    # --- Step 1: Calculate the Lorentz factor gamma_A ---
    # The Lorentz factor is defined as gamma = 1 / sqrt(1 - beta^2).
    try:
        gamma_A = 1 / math.sqrt(1 - beta_A**2)
    except ValueError:
        print("Error: beta_A must be less than 1.")
        return

    # --- Step 2: Calculate tan(theta) using the derived formula ---
    # The angle theta in the lab frame is given by the formula:
    # tan(theta) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))
    # This is derived from the Lorentz transformation of the momentum components of particle C.
    sqrt2 = math.sqrt(2)
    
    # We will output each number in the final equation.
    # The final equation is: tan_theta = 1 / (denominator)
    # where denominator = gamma_A * (1 + sqrt(2) * beta_A)
    
    term_in_parenthesis = 1 + sqrt2 * beta_A
    denominator = gamma_A * term_in_parenthesis
    tan_theta = 1 / denominator

    # --- Step 3: Calculate the angle theta in degrees ---
    # theta = arctan(tan_theta)
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)

    # --- Step 4: Print the results and the final answer ---
    print("Calculation of the angle of particle C in the lab frame.")
    print("-" * 50)
    print("Inputs and intermediate values for the equation:")
    print(f"  Velocity of parent particle A (beta_A) = {beta_A}")
    print(f"  Lorentz factor (gamma_A)              = {gamma_A}")
    print(f"  Value of sqrt(2)                     = {sqrt2}")
    print("\nFinal Equation for tan(theta):")
    print(f"  tan(theta) = 1 / ({gamma_A:.4f} * (1 + {sqrt2:.4f} * {beta_A}))")
    print(f"  tan(theta) = 1 / ({denominator:.4f})")
    print(f"  tan(theta) = {tan_theta}")
    print("-" * 50)
    print(f"\nThe angle in the lab frame is {theta_deg:.3f} degrees.")


if __name__ == "__main__":
    calculate_angle()
    # The final answer is now calculated and displayed above.
    # To conform to the output format, let's recalculate and format it here.
    beta_A = 0.95
    gamma_A = 1 / math.sqrt(1 - beta_A**2)
    tan_theta = 1 / (gamma_A * (1 + math.sqrt(2) * beta_A))
    theta_deg = math.degrees(math.atan(tan_theta))
    final_answer = round(theta_deg, 3)
    print(f'<<<{final_answer}>>>')