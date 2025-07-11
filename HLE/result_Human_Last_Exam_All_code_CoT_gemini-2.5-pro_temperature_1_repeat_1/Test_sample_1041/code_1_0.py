import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D simple random walk visits a target set
    before leaving a large disk, using a continuous potential theory approximation.
    """
    # Problem parameters
    R = 1000.0  # Radius of the disk
    z0 = (0, 300)  # Starting point
    a1 = (0, 0)   # First target point
    a2 = (2, 0)   # Second target point

    # Effective radius of a lattice point, derived from self-consistency
    epsilon = 0.25

    # --- Calculate terms for the probability formula ---

    # Distances from starting point to target points
    dist_z0_a1 = math.sqrt((z0[0] - a1[0])**2 + (z0[1] - a1[1])**2)
    dist_z0_a2 = math.sqrt((z0[0] - a2[0])**2 + (z0[1] - a2[1])**2)
    
    # Distance between target points
    dist_a1_a2 = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2)
    
    # The Green's function G_R(z, a) for z far from a and the boundary is approx log(R / |z-a|)
    # For z0, this approximation is excellent.
    num1 = math.log(R / dist_z0_a1)
    num2 = math.log(R / dist_z0_a2)
    numerator = num1 + num2
    
    # The denominator terms are derived from evaluating the solution on the target set
    # G_R(a1, a1) is singular, replaced by log(R/epsilon)
    den1 = math.log(R / epsilon)
    # G_R(a1, a2) is approx log(R / |a1-a2|)
    den2 = math.log(R / dist_a1_a2)
    denominator = den1 + den2
    
    # Calculate the final probability
    probability = numerator / denominator

    # --- Print the calculation step-by-step ---
    
    print("The probability P is calculated using the formula:")
    print("P = (log(R/|z0-a1|) + log(R/|z0-a2|)) / (log(R/epsilon) + log(R/|a1-a2|))")
    print("\nWhere:")
    print(f"  R = {R}")
    print(f"  z0 = {z0}, a1 = {a1}, a2 = {a2}")
    print(f"  |z0-a1| = {dist_z0_a1:.4f}")
    print(f"  |z0-a2| = {dist_z0_a2:.4f}")
    print(f"  |a1-a2| = {dist_a1_a2:.4f}")
    print(f"  epsilon = {epsilon}")

    print("\nNumerator calculation:")
    print(f"  log({R}/{dist_z0_a1:.4f}) = log({R/dist_z0_a1:.4f}) = {num1:.4f}")
    print(f"  log({R}/{dist_z0_a2:.4f}) = log({R/dist_z0_a2:.4f}) = {num2:.4f}")
    print(f"  Numerator = {num1:.4f} + {num2:.4f} = {numerator:.4f}")

    print("\nDenominator calculation:")
    print(f"  log({R}/{epsilon}) = log({R/epsilon:.4f}) = {den1:.4f}")
    print(f"  log({R}/{dist_a1_a2:.4f}) = log({R/dist_a1_a2:.4f}) = {den2:.4f}")
    print(f"  Denominator = {den1:.4f} + {den2:.4f} = {denominator:.4f}")

    print("\nFinal Probability:")
    print(f"  P = {numerator:.4f} / {denominator:.4f} = {probability:.4f}")
    
    # Format the answer to three significant digits
    print(f"\nThe probability to three significant digits is {probability:.3g}.")


solve_random_walk_probability()