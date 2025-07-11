import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem in ground effect.

    The problem is solved using the vortex mirror image method and thin aerofoil theory.
    The lift ratio L1/L2 is found to be (1-A)/(1+A), where A is an interaction factor
    dependent on the geometry of the system.
    """
    # We can set c=1 without loss of generality as the result is a dimensionless ratio.
    c = 1.0
    s = 0.5 * c  # Separation distance
    h = 0.5 * c  # Ride height

    # The interaction factor 'A' is derived from the induced velocity calculations.
    # A = (c/2) * [1/s - s / (s^2 + 4*h^2)]
    
    term1 = 1 / s
    term2 = s / (s**2 + 4 * h**2)
    
    A = (c / 2) * (term1 - term2)

    # The lift ratio L1/L2 is given by (1-A)/(1+A)
    # assuming the geometric angles of attack of the two aerofoils are the same.
    if (1 + A) == 0:
        # Avoid division by zero, though it's not expected for this problem
        lift_ratio = float('inf')
    else:
        lift_ratio = (1 - A) / (1 + A)

    # Print the results, showing the numbers in the final equation.
    print(f"Given geometry:")
    print(f"  Chord c")
    print(f"  Separation s = 0.5 * c")
    print(f"  Ride height h = 0.5 * c\n")
    
    print(f"The interaction factor 'A' is calculated as:")
    print(f"A = (c/2) * [1/s - s / (s^2 + 4*h^2)]")
    print(f"A = ({c}/2) * [1/{s} - {s} / ({s}^2 + 4*{h}^2)] = {A:.4f}\n")

    print(f"The lift ratio L1/L2 is calculated using the formula (1 - A) / (1 + A):")
    numerator = 1 - A
    denominator = 1 + A
    print(f"L1/L2 = (1 - {A:.4f}) / (1 + {A:.4f})")
    print(f"L1/L2 = {numerator:.4f} / {denominator:.4f}")
    print(f"L1/L2 = {lift_ratio:.4f}\n")
    
    # Also show the exact fractional answer
    # A = 4/5, so ratio = (1 - 4/5) / (1 + 4/5) = (1/5) / (9/5) = 1/9
    print(f"The exact fractional result is 1/9.")

calculate_lift_ratio()
<<<0.1111>>>