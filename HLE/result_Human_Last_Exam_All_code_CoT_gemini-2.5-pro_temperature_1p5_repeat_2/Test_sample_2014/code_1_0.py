def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem in ground effect.

    The problem is solved using the mirror image method from potential flow theory.
    The final ratio depends on a geometric interaction factor 'K'.
    """

    # We can set the chord length 'c' to 1.0 since it is a reference length
    # and will cancel out in the calculation of the non-dimensional factor K.
    c = 1.0

    # Define the separation distance 's' and ride height 'h' as given in the problem.
    s = 0.5 * c
    h = 0.5 * c

    print("--- Step 1: Calculate the interaction factor K ---")
    print(f"The problem specifies a chord 'c', a separation s = c/2, and a ride height h = c/2.")
    print(f"Using a reference chord c = {c}, we have s = {s} and h = {h}.")
    print("\nThe interaction factor K is derived from the induced velocities and is given by the formula:")
    print("K = (c / (2*s)) * (4*h^2) / (s^2 + 4*h^2)")

    # Calculate the individual components of the formula for K to show the steps.
    term1 = c / (2 * s)
    term2_num = 4 * h**2
    term2_den = s**2 + 4 * h**2
    term2 = term2_num / term2_den
    K = term1 * term2

    print("\nPlugging in the values:")
    print(f"K = ({c} / (2 * {s})) * (4 * {h:.1f}^2) / ({s:.1f}^2 + 4 * {h:.1f}^2)")
    print(f"K = ({term1}) * ({term2_num}) / ({s**2} + {4 * h**2})")
    print(f"K = ({term1}) * ({term2_num}) / ({term2_den})")
    print(f"K = {term1} * {term2:.2f}")
    print(f"So, K = {K:.2f}")

    print("\n--- Step 2: Calculate the final lift ratio L1/L2 ---")
    print("The lift ratio L1/L2 is related to K by the formula:")
    print("L1/L2 = (1 + K) / (1 - K)")

    # Calculate the numerator and denominator for the final ratio.
    numerator = 1 + K
    denominator = 1 - K
    lift_ratio = numerator / denominator

    print("\nPlugging in the value of K:")
    print(f"L1/L2 = (1 + {K:.2f}) / (1 - {K:.2f})")
    print(f"L1/L2 = {numerator:.2f} / {denominator:.2f}")
    print(f"\nFinal Lift Ratio L1/L2 = {lift_ratio:.1f}")

# Execute the calculation function.
calculate_lift_ratio()