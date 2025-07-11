import numpy as np
import argparse

def calculate_l(d, lmbda):
    """
    Calculates the value of l(d, lambda) based on the derived formula.

    The core of the problem is to find the logarithm of a ratio of probability densities.
    After analyzing the provided sampling function and assuming a standard geometric
    interpretation (the exponential map on a sphere), the density p(x) for a point x
    on the sphere is found to be proportional to exp(- (1/(2*lambda)) * ||v||^2),
    where v is the representation of x in the tangent space at a specific pole p0.
    The squared norm of this vector v is given by ||v||^2 = arccos(x^T * p0)^2.

    This leads to the final formula:
    l(d, lambda) = (1 / (2 * lambda)) * (arccos(x2^T * p0)^2 - arccos(x1^T * p0)^2)
    """
    # Validate inputs based on the problem definition
    if not (isinstance(d, int) and d >= 4):
        raise ValueError("d must be an integer greater than or equal to 4.")
    if not (isinstance(lmbda, (int, float)) and lmbda >= 1):
        raise ValueError("lambda must be a number greater than or equal to 1.")

    # Define the pole p0 = (1/sqrt(d)) * [1, 1, ..., 1]
    p0 = np.ones(d) / np.sqrt(d)

    # Define the vector x1 = (e1 + e2 + e3) / sqrt(3)
    x1 = np.zeros(d)
    x1[:3] = 1.0
    x1 /= np.sqrt(3.0)

    # Define the vector x2 = (e3 + e4) / sqrt(2)
    x2 = np.zeros(d)
    # Python uses 0-based indexing, so e3 and e4 are at indices 2 and 3
    x2[2:4] = 1.0
    x2 /= np.sqrt(2.0)

    # Calculate the dot products of xi and p0
    dot1 = np.dot(x1, p0)
    dot2 = np.dot(x2, p0)
    
    # Clip results to the valid domain for arccos, [-1, 1], to prevent
    # floating-point errors from causing a domain error.
    dot1_clipped = np.clip(dot1, -1.0, 1.0)
    dot2_clipped = np.clip(dot2, -1.0, 1.0)

    # Calculate theta = arccos(dot_product) for each vector
    theta1 = np.arccos(dot1_clipped)
    theta2 = np.arccos(dot2_clipped)
    
    # The squared norm of v is theta^2
    v1_norm_sq = theta1**2
    v2_norm_sq = theta2**2

    # Apply the final formula for l(d, lambda)
    result = (1 / (2 * lmbda)) * (v2_norm_sq - v1_norm_sq)

    # Output the detailed calculation steps as requested
    print("--- Calculation Details ---")
    print(f"l(d, lambda) = (1 / (2 * lambda)) * (arccos(x2^T * p0)^2 - arccos(x1^T * p0)^2)")
    print(f"\nFor d={d}, lambda={lmbda}:")
    
    # Step 1: Dot products
    print(f"1. Dot products:")
    print(f"   x1^T * p0 = {dot1:.8f}")
    print(f"   x2^T * p0 = {dot2:.8f}")
    
    # Step 2: Squared norms of v (theta^2)
    print(f"\n2. Squared norms of tangent vectors (||v||^2 = arccos^2):")
    print(f"   ||v1||^2 = arccos({dot1_clipped:.8f})^2 = {theta1:.8f}^2 = {v1_norm_sq:.8f}")
    print(f"   ||v2||^2 = arccos({dot2_clipped:.8f})^2 = {theta2:.8f}^2 = {v2_norm_sq:.8f}")
    
    # Step 3: Final computation
    print(f"\n3. Final result for l({d}, {lmbda}):")
    term_a = 1 / (2 * lmbda)
    term_b = v2_norm_sq - v1_norm_sq
    print(f"   l({d}, {lmbda}) = ({term_a:.8f}) * ({v2_norm_sq:.8f} - {v1_norm_sq:.8f})")
    print(f"   l({d}, {lmbda}) = ({term_a:.8f}) * ({term_b:.8f})")
    print(f"   l({d}, {lmbda}) = {result:.8f}")
    print("-------------------------\n")

    # Final answer in the requested format
    print("<<<" + f"{result}" + ">>>")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate l(d, lambda) based on the provided sampling function.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--d',
        type=int,
        default=4,
        help='Dimension d, must be an integer >= 4. Default: 4'
    )
    parser.add_argument(
        '--lmbda',
        type=float,
        default=1.0,
        help='Lambda parameter, must be a float >= 1.0. Default: 1.0'
    )
    args = parser.parse_args()
    
    try:
        calculate_l(args.d, args.lmbda)
    except ValueError as e:
        print(f"Error: {e}")
