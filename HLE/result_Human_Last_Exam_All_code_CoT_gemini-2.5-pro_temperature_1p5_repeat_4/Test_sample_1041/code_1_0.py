import numpy as np

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D simple random walk visits a target set
    before leaving a large disk, using a continuous potential theory approximation.
    """
    # Parameters of the problem
    R = 1000.0  # Radius of the disk
    start_pos = np.array([0, 300])  # Starting position of the random walk
    target_pos_1 = np.array([0, 0])   # First point in the target set
    target_pos_2 = np.array([2, 0])   # Second point in the target set

    # In the continuous approximation, we replace the discrete target points with small disks.
    # A natural choice for the effective radius of a lattice point is the lattice spacing itself.
    r0 = 1.0

    # Calculate geometric quantities
    d = np.linalg.norm(target_pos_1 - target_pos_2)
    target_center = (target_pos_1 + target_pos_2) / 2.0
    rs = np.linalg.norm(start_pos - target_center)

    # The hitting probability can be estimated using a formula derived from potential theory.
    # P = [2 * ln(R/rs)] / [ln(R/r0) + ln(R/d)]
    # This formula represents the electrostatic potential at the starting point 'start_pos'
    # in a system where two small disks (representing the target points) are held
    # at potential 1, and a large circle (representing the disk boundary) is grounded (potential 0).
    
    numerator = 2 * np.log(R / rs)
    denominator = np.log(R / r0) + np.log(R / d)
    
    probability = numerator / denominator

    # Output the parameters and the final result
    print("Problem Parameters and Calculation:")
    print(f"Radius of the disk, R = {R}")
    print(f"Distance from start to target center, rs = {rs:.4f}")
    print(f"Effective radius of a target point, r0 = {r0}")
    print(f"Distance between target points, d = {d}")
    
    print("\nThe formula used is P = [2 * ln(R/rs)] / [ln(R/r0) + ln(R/d)]")
    
    # Printing the equation with the numerical values
    print("\nPlugging in the numbers:")
    print(f"P = [2 * ln({R:.1f}/{rs:.4f})] / [ln({R:.1f}/{r0:.1f}) + ln({R:.1f}/{d:.1f})]")
    print(f"P = [2 * ln({R/rs:.4f})] / [ln({R/r0:.4f}) + ln({R/d:.4f})]")
    print(f"P = [{numerator:.4f}] / [{np.log(R/r0):.4f} + {np.log(R/d):.4f}]")
    print(f"P = {numerator:.4f} / {denominator:.4f}")

    print(f"\nThe calculated probability is: {probability:.5f}")
    print(f"The probability rounded to three significant digits is: {probability:.3g}")


solve_random_walk_probability()
<<<0.183>>>