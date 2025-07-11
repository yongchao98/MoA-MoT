import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.

    This function follows these steps:
    1.  Sets a nominal chord length, c=1.0, as the ratio is dimensionless.
    2.  Defines the separation 's' and height 'h' based on the problem statement.
    3.  Calculates the influence coefficient 'K' using the derived formula:
        K = (2 * c * h^2) / (s * (s^2 + 4*h^2))
    4.  Calculates the final lift ratio L1/L2 using the formula:
        L1/L2 = (1 + K) / (1 - K)
    5.  Prints each step of the calculation, showing the intermediate values.
    """
    # 1. Define chord length and parameters
    c = 1.0
    s = 0.5 * c
    h = 0.5 * c

    print("--- Problem Setup ---")
    print(f"Aerofoil chord length, c = {c}")
    print(f"Separation distance, s = 1/2c = {s}")
    print(f"Ride height, h = 1/2c = {h}")
    print("-" * 25)

    # 2. Calculate the influence coefficient K
    print("--- Calculating Influence Coefficient (K) ---")
    # Formula: K = (2 * c * h^2) / (s * (s^2 + 4*h^2))
    k_numerator = 2 * c * h**2
    k_denominator = s * (s**2 + 4 * h**2)
    K = k_numerator / k_denominator

    print(f"Formula for K = (2 * c * h^2) / (s * (s^2 + 4*h^2))")
    print(f"Numerator = 2 * {c:.1f} * {h:.1f}^2 = {k_numerator:.3f}")
    print(f"Denominator = {s:.1f} * ({s:.1f}^2 + 4 * {h:.1f}^2) = {k_denominator:.3f}")
    print(f"K = {k_numerator:.3f} / {k_denominator:.3f} = {K:.2f}")
    print("-" * 25)

    # 3. Calculate the lift ratio L1/L2
    print("--- Calculating Lift Ratio (L1/L2) ---")
    # Formula: L1/L2 = (1 + K) / (1 - K)
    ratio_numerator = 1 + K
    ratio_denominator = 1 - K
    lift_ratio = ratio_numerator / ratio_denominator

    print(f"Formula for L1/L2 = (1 + K) / (1 - K)")
    # Final equation with numbers
    print(f"L1/L2 = (1 + {K:.2f}) / (1 - {K:.2f})")
    print(f"L1/L2 = {ratio_numerator:.2f} / {ratio_denominator:.2f}")
    print(f"Final Lift Ratio L1/L2 = {lift_ratio:.1f}")


if __name__ == "__main__":
    calculate_lift_ratio()
    final_answer = 9.0 # Calculated value
    # The final answer will be printed in the required format by the execution environment.