import numpy as np

def solve_task():
    """
    Calculates the value of l(d, lambda) based on the derived formula and prints
    a step-by-step breakdown of the calculation.

    The user can configure the values for d and lambda below.
    """
    # --- User-configurable values ---
    # d: dimension, must be an integer >= 4
    # lambd: lambda parameter, must be >= 1
    d = 4
    lambd = 1.0
    # ------------------------------------

    # Validate inputs
    if not isinstance(d, int) or d < 4:
        print(f"Error: d must be an integer greater than or equal to 4. Current value is {d}.")
        return
    if not isinstance(lambd, (int, float)) or lambd < 1:
        print(f"Error: lambda must be a number greater than or equal to 1. Current value is {lambd}.")
        return

    # Derived formula for l(d, lambda):
    # l = ln[g(v1)/g(v2)] where g(v) = exp(-||v||^2/(2*lambda)) * (sin(||v||)/||v||)^(d-2)
    # v1 and v2 are the principal pre-images of x1 and x2.

    # Step 1: Calculate the norms of the pre-images, theta_1 and theta_2
    theta_1 = np.arccos(np.sqrt(3 / d))
    theta_2 = np.arccos(np.sqrt(2 / d))

    # Step 2: Calculate the two main components of the log-likelihood expression
    # Term 1: Comes from the exponential part of g(v)
    term1 = (theta_2**2 - theta_1**2) / (2 * lambd)

    # Term 2: Comes from the sinc part of g(v)
    # sinc(x) = sin(x)/x
    sinc_theta1 = np.sin(theta_1) / theta_1
    sinc_theta2 = np.sin(theta_2) / theta_2
    
    # The argument to the logarithm should be positive. This holds for d>=4.
    term2 = (d - 2) * np.log(sinc_theta1 / sinc_theta2)
    
    # Final result is the sum of the two terms
    result = term1 + term2

    # Print the detailed breakdown of the calculation
    print("### Calculation of l(d, lambda) ###\n")
    print(f"This script calculates l(d, lambda) for d = {d} and lambda = {lambd}.\n")
    print("The final formula is: l(d, lambda) = Term1 + Term2, where")
    print("Term1 = (theta_2^2 - theta_1^2) / (2 * lambda)")
    print("Term2 = (d - 2) * ln[ (sin(theta_1)/theta_1) / (sin(theta_2)/theta_2) ]\n")
    print("-" * 50)
    
    print(f"Step 1: Calculate theta_1 and theta_2")
    print(f"theta_1 = arccos(sqrt(3 / {d})) = {theta_1:.6f}")
    print(f"theta_2 = arccos(sqrt(2 / {d})) = {theta_2:.6f}\n")

    print(f"Step 2: Calculate Term 1")
    print(f"Term1 = (({theta_2:.6f})^2 - ({theta_1:.6f})^2) / (2 * {lambd})")
    print(f"      = {term1:.6f}\n")

    print(f"Step 3: Calculate Term 2")
    print(f"sin(theta_1)/theta_1 = {sinc_theta1:.6f}")
    print(f"sin(theta_2)/theta_2 = {sinc_theta2:.6f}")
    print(f"Term2 = ({d} - 2) * ln({sinc_theta1:.6f} / {sinc_theta2:.6f})")
    print(f"      = {d - 2} * ln({sinc_theta1 / sinc_theta2:.6f})")
    print(f"      = {term2:.6f}\n")
    
    print("-" * 50)
    print(f"Final Result: l({d}, {lambd}) = Term1 + Term2")
    print(f"l({d}, {lambd}) = {term1:.6f} + {term2:.6f} = {result:.6f}")

# Execute the calculation
solve_task()