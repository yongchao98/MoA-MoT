import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D random walk visits a target set before
    leaving a large disk, using a continuous potential approximation.
    """
    # Parameters from the problem statement
    z0 = (0, 300)
    A = [(0, 0), (2, 0)]
    R = 1000

    # 1. Calculate the center of the target set A
    zc = ((A[0][0] + A[1][0]) / 2, (A[0][1] + A[1][1]) / 2)

    # 2. Calculate the distance d between the points in A
    d = math.sqrt((A[1][0] - A[0][0])**2 + (A[1][1] - A[0][1])**2)

    # 3. Calculate the effective radius of the target set A
    # The effective radius of a single lattice point, r0, is taken as 1.
    r0 = 1
    r_eff = math.sqrt(r0 * d)

    # 4. Calculate the distance from the starting point to the center of A
    dist_z0_zc = math.sqrt((z0[0] - zc[0])**2 + (z0[1] - zc[1])**2)

    # 5. Calculate the probability using the formula
    log_R = math.log(R)
    log_dist = math.log(dist_z0_zc)
    log_reff = math.log(r_eff)

    numerator = log_R - log_dist
    denominator = log_R - log_reff
    probability = numerator / denominator

    # 6. Print the equation with the numbers and the result
    print("The formula for the probability is P = (ln(R) - ln(|z - zc|)) / (ln(R) - ln(r_eff))")
    print("\nSubstituting the given values:")
    print(f"R = {R}")
    print(f"z = {z0}")
    print(f"A = {A}")
    print(f"zc = {zc}")
    print(f"|z - zc| = {dist_z0_zc:.4f}")
    print(f"d = {d:.4f}")
    print(f"r_eff = {r_eff:.4f}")
    
    print("\nPlugging numbers into the formula:")
    print(f"P = (ln({R}) - ln({dist_z0_zc:.4f})) / (ln({R}) - ln({r_eff:.4f}))")
    print(f"P = ({log_R:.4f} - {log_dist:.4f}) / ({log_R:.4f} - {log_reff:.4f})")
    print(f"P = {numerator:.4f} / {denominator:.4f}")
    print(f"\nThe calculated probability is: {probability:.5f}")
    print(f"The probability rounded to three significant digits is: {probability:.3f}")


solve_random_walk_probability()