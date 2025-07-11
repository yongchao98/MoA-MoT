import numpy as np

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D random walk starting at (0, 300)
    visits the set {(0,0), (2,0)} before leaving the disk of radius 1000.
    """
    # Parameters
    R = 1000.0
    epsilon = 1.0  # Lattice spacing, effective radius of target points
    
    # Complex number representation of the points
    z0 = 0 + 300j # Start point (0, 300)
    z1 = 0 + 0j   # Target point 1 (0, 0)
    z2 = 2 + 0j   # Target point 2 (2, 0)

    # Calculate distances between points
    d01 = np.abs(z0 - z1)
    d02 = np.abs(z0 - z2)
    d12 = np.abs(z1 - z2)

    # Numerator terms of the formula
    num1 = np.log(R / d01)
    num2 = np.log(R / d02)
    numerator = num1 + num2
    
    # Denominator terms of the formula
    den1 = np.log(R / epsilon)
    den2 = np.log(R / d12)
    denominator = den1 + den2

    # Calculate the probability
    probability = numerator / denominator

    # Output the equation with the calculated values
    print("The probability p is calculated using the formula:")
    print("p = (log(R/|z0-z1|) + log(R/|z0-z2|)) / (log(R/epsilon) + log(R/|z1-z2|))")
    print("\nPlugging in the values:")
    print(f"R = {R}")
    print(f"z0 = {z0}, z1 = {z1}, z2 = {z2}")
    print(f"epsilon = {epsilon}")
    print(f"|z0-z1| = |{z0 - z1}| = {d01:.4f}")
    print(f"|z0-z2| = |{z0 - z2}| = {d02:.4f}")
    print(f"|z1-z2| = |{z1 - z2}| = {d12:.4f}")
    print("\nEquation with numbers:")
    print(f"p = (log({R}/{d01:.4f}) + log({R}/{d02:.4f})) / (log({R}/{epsilon:.4f}) + log({R}/{d12:.4f}))")
    print(f"p = ({num1:.4f} + {num2:.4f}) / ({den1:.4f} + {den2:.4f})")
    print(f"p = {numerator:.4f} / {denominator:.4f}")
    
    # Final result
    print(f"\nThe probability is: {probability:.3f}")

solve_random_walk_probability()