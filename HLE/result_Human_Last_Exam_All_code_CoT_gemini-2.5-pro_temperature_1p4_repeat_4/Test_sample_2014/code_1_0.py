import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    and in ground effect using the mirror image method.
    """
    # For calculation purposes, we can normalize the chord length c to 1,
    # as it will cancel out in the dimensionless factor 'A'.
    c = 1.0

    # Given geometric parameters
    s = 0.5 * c  # Separation between aerofoils
    h = 0.5 * c  # Ride height above the ground

    print("Step 1: Define problem parameters (normalized with chord c=1)")
    print(f"Chord, c = {c}")
    print(f"Separation, s = {s}")
    print(f"Ride Height, h = {h}")
    print("-" * 30)

    # Step 2: Calculate the dimensionless interaction factor 'A'.
    # The lift equations are:
    # Γ1 = Γ0 - A * Γ2
    # Γ2 = Γ0 + A * Γ1
    # where A = (c/2) * [1/s - s/(s^2 + 4h^2)]
    
    bracket_term_1 = 1 / s
    bracket_term_2 = s / (s**2 + 4 * h**2)
    A = (c / 2) * (bracket_term_1 - bracket_term_2)

    print("Step 2: Calculate the interaction factor 'A'")
    print(f"A = (c/2) * [1/s - s/(s^2 + 4*h^2)]")
    print(f"A = ({c}/2) * [1/{s} - {s}/({s}^2 + 4*{h}^2)]")
    print(f"A = {A:.4f}")
    print("-" * 30)

    # Step 3: Solve for the lift ratio L1/L2.
    # The solution to the system of equations gives:
    # L1/L2 = Γ1/Γ2 = (1 - A) / (1 + A)
    
    numerator = 1 - A
    denominator = 1 + A
    lift_ratio = numerator / denominator

    print("Step 3: Calculate the final lift ratio L1/L2")
    print(f"The formula for the lift ratio is: L1/L2 = (1 - A) / (1 + A)")
    # Using f-string to format the output to show the numbers in the equation
    print("Plugging in the value for A:")
    print(f"L1/L2 = (1 - {A:.2f}) / (1 + {A:.2f}) = {numerator:.2f} / {denominator:.2f}")
    print(f"Final calculated ratio: {lift_ratio:.4f}")
    print(f"The exact fractional result is 1/9.")

calculate_lift_ratio()