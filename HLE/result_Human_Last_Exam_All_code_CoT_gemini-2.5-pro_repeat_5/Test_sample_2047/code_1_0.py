import numpy as np
import argparse

def calculate_l(d, lambd):
    """
    Calculates the value of l(d, lambda) based on the derived formula.

    The formula is: l(d, λ) = (1 / (2λ)) * [ (arccos(sqrt(2/d)))^2 - (arccos(sqrt(3/d)))^2 ]

    Args:
        d (int): The dimension, must be >= 4.
        lambd (float): The lambda parameter, must be >= 1.

    Returns:
        float: The calculated value of l(d, lambda).
    """
    if not (isinstance(d, int) and d >= 4):
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not (isinstance(lambd, (int, float)) and lambd >= 1):
        raise ValueError("lambda must be a number greater than or equal to 1.")

    # Calculate the terms inside the arccos function
    val_for_theta1 = np.sqrt(3.0 / d)
    val_for_theta2 = np.sqrt(2.0 / d)

    # Calculate the squared angles
    theta1_sq = np.arccos(val_for_theta1)**2
    theta2_sq = np.arccos(val_for_theta2)**2
    
    # Calculate the coefficient
    coeff = 1.0 / (2.0 * lambd)

    # Calculate the final value of l(d, lambda)
    l_value = coeff * (theta2_sq - theta1_sq)

    # Print the equation with the computed numbers, as requested.
    print(f"Calculating l(d, λ) for d={d}, λ={lambd}:")
    print(f"l({d}, {lambd}) = (1 / (2 * {lambd})) * [ (arccos(sqrt(2/{d})))^2 - (arccos(sqrt(3/{d})))^2 ]")
    print(f"l({d}, {lambd}) = ({coeff:.4f}) * [ (arccos({val_for_theta2:.4f}))^2 - (arccos({val_for_theta1:.4f}))^2 ]")
    print(f"l({d}, {lambd}) = ({coeff:.4f}) * [ {theta2_sq:.4f} - {theta1_sq:.4f} ]")
    print(f"Final Result: l({d}, {lambd}) = {l_value}")
    
    return l_value

if __name__ == '__main__':
    # Example usage:
    # You can run this script from the command line with arguments, e.g.:
    # python your_script_name.py --d 5 --lambd 2
    # If no arguments are provided, it uses the default example values.
    parser = argparse.ArgumentParser(description="Calculate l(d, lambda).")
    parser.add_argument('--d', type=int, default=5, help='Dimension d (>= 4)')
    parser.add_argument('--lambd', type=float, default=2.0, help='Lambda parameter (>= 1)')
    
    args = parser.parse_args()
    
    try:
        calculate_l(args.d, args.lambd)
    except ValueError as e:
        print(f"Error: {e}")
