import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.
    """
    # For calculation purposes, we can set the chord length 'c' to 1.0,
    # as the final ratio is independent of its value.
    c = 1.0
    print(f"Assuming a chord length, c = {c}\n")

    # Define the separation 's' and ride height 'h' based on the problem statement.
    s = 0.5 * c
    h = 0.5 * c
    print(f"Given separation s = 1/2*c = {s}")
    print(f"Given ride height h = c/2 = {h}\n")

    # The lift ratio L1/L2 can be found using the formula derived from the mirror image method:
    # L1/L2 = (1/c + K) / (1/c - K), where K = 2*h^2 / (s*(s^2 + 4*h^2))
    print("The formula for the lift ratio is L1/L2 = (1/c + K) / (1/c - K)")
    print("where the interaction term K = (2*h^2) / (s*(s^2 + 4*h^2))\n")

    # Step 1: Calculate the interaction term K
    term_k_numerator = 2 * h**2
    term_k_denominator = s * (s**2 + 4 * h**2)
    k = term_k_numerator / term_k_denominator
    
    print("First, calculating the value of K:")
    print(f"K = (2 * {h}^2) / ({s} * ({s}^2 + 4 * {h}^2))")
    print(f"K = ({term_k_numerator}) / ({term_k_denominator})")
    print(f"K = {k:.4f}\n")

    # Step 2: Substitute K into the lift ratio formula
    ratio_numerator = (1 / c) + k
    ratio_denominator = (1 / c) - k

    print("Next, substituting K into the lift ratio equation:")
    print(f"L1/L2 = ({1/c:.4f} + {k:.4f}) / ({1/c:.4f} - {k:.4f})")
    print(f"L1/L2 = {ratio_numerator:.4f} / {ratio_denominator:.4f}\n")

    # Step 3: Calculate the final result
    lift_ratio = ratio_numerator / ratio_denominator
    
    print(f"The final calculated lift ratio L1/L2 is: {lift_ratio}")

# Run the calculation
calculate_lift_ratio()