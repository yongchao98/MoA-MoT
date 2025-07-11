import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    in ground effect using the mirror image method.
    """
    # We can set the chord length c to 1.0, as the result is a dimensionless ratio.
    c = 1.0

    # Define the separation distance s and ride height h based on the problem statement.
    s = 0.5 * c
    h = 0.5 * c

    print("--- Problem Parameters ---")
    print(f"Chord length, c = {c}")
    print(f"Separation distance, s = 1/2 * c = {s}")
    print(f"Ride height, h = 1/2 * c = {h}\n")

    # The interaction between the aerofoils and their ground images is captured by a
    # dimensionless influence coefficient, M.
    # The formula for M is: M = (c/2) * [1/s - s / (s^2 + 4*h^2)]
    
    print("--- Calculation of Influence Coefficient (M) ---")
    
    # Calculate M step-by-step
    s_squared = s**2
    h_squared = h**2
    denominator = s_squared + 4 * h_squared
    term1 = 1 / s
    term2 = s / denominator
    M = (c / 2) * (term1 - term2)
    
    print(f"Equation for M: M = (c/2) * [1/s - s / (s^2 + 4*h^2)]")
    print(f"Substituting values: M = ({c}/2) * [1/{s} - {s} / ({s}^2 + 4*{h}^2)]")
    print(f"Calculating denominator: s^2 + 4*h^2 = {s_squared} + 4*{h_squared} = {denominator}")
    print(f"Calculating M: M = {c/2} * [{term1} - {s}/{denominator}]")
    print(f"M = {c/2} * [{term1} - {term2}]")
    print(f"M = {M:.4f}\n")
    
    # The lift ratio L1/L2 can be calculated from M using the formula: L1/L2 = (1 - M) / (1 + M)
    print("--- Calculation of Lift Ratio (L1/L2) ---")
    numerator_lift_ratio = 1 - M
    denominator_lift_ratio = 1 + M
    lift_ratio = numerator_lift_ratio / denominator_lift_ratio
    
    print(f"Equation for lift ratio: L1/L2 = (1 - M) / (1 + M)")
    print(f"Substituting the value of M: L1/L2 = (1 - {M:.4f}) / (1 + {M:.4f})")
    print(f"Final Equation: L1/L2 = {numerator_lift_ratio:.4f} / {denominator_lift_ratio:.4f}")
    
    print("\n--- Result ---")
    print(f"The final calculated lift ratio L1/L2 is: {lift_ratio:.4f}")
    
    # Representing as a simple fraction
    # The result 0.1111 is 1/9
    print(f"As a simplified fraction, the result is 1/9.")

# Execute the calculation
calculate_lift_ratio()