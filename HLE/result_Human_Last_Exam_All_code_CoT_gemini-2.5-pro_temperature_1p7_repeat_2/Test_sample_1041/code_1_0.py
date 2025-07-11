import math

def solve_random_walk_probability():
    """
    Calculates the probability for the given random walk problem
    using an approximation from potential theory.
    """
    # Define the parameters from the problem description
    R = 1000.0  # Radius of the boundary disk
    start_pos = (0, 300)
    target_point1 = (0, 0)
    target_point2 = (2, 0)

    # The center of the target set is the midpoint of the two target points
    target_center = (
        (target_point1[0] + target_point2[0]) / 2.0,
        (target_point1[1] + target_point2[1]) / 2.0
    )

    # The distance 'd' between the two points in the target set
    d = math.sqrt((target_point1[0] - target_point2[0])**2 + (target_point1[1] - target_point2[1])**2)
    
    # The effective radius 'r_A' of the target set, approximated by the capacitary radius d/2
    r_A = d / 2.0

    # The distance 'd_S' from the starting point to the center of the target set
    d_S = math.sqrt((start_pos[0] - target_center[0])**2 + (start_pos[1] - target_center[1])**2)

    # Apply the potential theory formula for the probability P
    # P = log(R / d_S) / log(R / r_A)
    numerator = math.log(R / d_S)
    denominator = math.log(R / r_A)
    probability = numerator / denominator

    # Output the steps of the calculation
    print("Calculating the probability using the formula: P = log(R / d_S) / log(R / r_A)")
    print(f"R (boundary radius) = {R}")
    print(f"d_S (distance from start to target center) = {d_S:.3f}")
    print(f"r_A (effective target radius) = {r_A}")
    print("\nFinal equation with numbers:")
    print(f"P = log({R} / {d_S:.3f}) / log({R} / {r_A})")
    print(f"P = log({R/d_S:.4f}) / log({R/r_A:.1f})")
    print(f"P = {numerator:.4f} / {denominator:.4f}")
    
    # Final answer
    print(f"\nThe calculated probability is approximately: {probability:.4f}")
    print(f"The answer with three significant digits is {probability:.3g}")

solve_random_walk_probability()