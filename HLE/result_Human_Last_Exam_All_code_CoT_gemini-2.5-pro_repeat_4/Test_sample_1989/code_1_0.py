import math

def calculate_corrector(r, theta, A, B):
    """
    Calculates the corrector term for the large-distance behavior of the solution.

    The corrector C is given by the formula: C = r^P * e^E
    where P = A*(1-cos(theta)) + B*sin(theta) and E = B*theta.

    Args:
        r (float): The distance from the origin (must be > 0).
        theta (float): The angle in radians.
        A (float): The coefficient A from the problem description.
        B (float): The coefficient B from the problem description.

    Returns:
        float: The value of the corrector term.
    """
    if r <= 0:
        raise ValueError("r must be positive.")

    # Calculate the exponent of r
    power_of_r = A * (1 - math.cos(theta)) + B * math.sin(theta)
    
    # Calculate the exponent of e
    exp_argument = B * theta
    
    # Calculate the r-term and the exponential term
    r_term = r ** power_of_r
    exp_term = math.exp(exp_argument)

    # The final corrector is the product of the two
    corrector = r_term * exp_term

    print("--- Calculation Details ---")
    print(f"Parameters: r={r}, theta={theta:.4f}, A={A}, B={B}")
    print(f"Base 'r' of the power term: {r}")
    print(f"Exponent 'P' of the r-term: A(1-cos(theta)) + B*sin(theta) = {power_of_r:.4f}")
    print(f"Base 'e' of the exponential term: {math.e:.4f}")
    print(f"Exponent 'E' of the e-term: B*theta = {exp_argument:.4f}")
    print(f"Final Corrector Value = r^P * e^E = {corrector:.4f}")
    
    return corrector

if __name__ == '__main__':
    # Example usage of the function
    # Let's choose some example values for the parameters
    r_val = 10.0
    theta_val = math.pi / 2  # 90 degrees
    A_val = 1.0
    B_val = 1.0

    print("Calculating the corrector for a sample case...")
    final_corrector = calculate_corrector(r_val, theta_val, A_val, B_val)
    print("\n--- Final Result ---")
    # The final equation for the corrector is C = r^(A(1-cos(theta)) + B*sin(theta)) * e^(B*theta)
    # Let's print the numbers in this equation for our example
    print(f"The corrector for r={r_val}, theta=pi/2, A={A_val}, B={B_val} is: {final_corrector:.4f}")
