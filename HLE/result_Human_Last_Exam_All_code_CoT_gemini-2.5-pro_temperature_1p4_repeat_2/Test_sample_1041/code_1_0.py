import math

def solve_random_walk_problem():
    """
    Calculates the probability that a 2D random walk visits a target set
    before leaving a large disk.
    """
    # Problem parameters
    z0 = (0, 300)
    a1 = (0, 0)
    a2 = (2, 0)
    R = 1000.0  # Radius of the disk
    a = 1.0     # Effective radius of a single lattice point

    print("Problem Parameters:")
    print(f"Start point z0: {z0}")
    print(f"Target set A: {{{a1}, {a2}}}")
    print(f"Disk radius R: {R}")
    print("-" * 30)

    # Step 1: Calculate relevant distances
    dist_z0_a1 = math.sqrt((z0[0] - a1[0])**2 + (z0[1] - a1[1])**2)
    dist_z0_a2 = math.sqrt((z0[0] - a2[0])**2 + (z0[1] - a2[1])**2)
    dist_a1_a2 = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2)

    # Step 2: Define the single-point hitting probability function components
    # p(z, w) is proportional to log(R / |z-w|)
    # The normalization factor is log(R/a)
    log_R_a = math.log(R / a)

    # p_ij = p(a_i, a_j) used for finding the coefficients C_i
    # p(a_i, a_i) is 1 by definition.
    # p(a_1, a_2) = p(a_2, a_1) due to symmetry.
    p_a1_a2 = math.log(R / dist_a1_a2) / log_R_a
    
    # Step 3: Solve for coefficients C1 and C2
    # The system of equations is:
    # C1 + C2 * p_a1_a2 = 1
    # C1 * p_a1_a2 + C2 = 1
    # By symmetry, C1 = C2.
    # C1 * (1 + p_a1_a2) = 1  => C1 = 1 / (1 + p_a1_a2)
    C = 1 / (1 + p_a1_a2)
    C1 = C
    C2 = C
    
    print("Calculating coefficients C1 and C2:")
    print(f"Interaction term p(a1, a2) = log({R}/{dist_a1_a2}) / log({R}/{a}) = {p_a1_a2:.4f}")
    print(f"Solving C*(1 + {p_a1_a2:.4f}) = 1  => C1 = C2 = {C:.4f}")
    print("-" * 30)

    # Step 4: Calculate the hitting probabilities from z0 to each point in A
    p_z0_a1 = math.log(R / dist_z0_a1) / log_R_a
    p_z0_a2 = math.log(R / dist_z0_a2) / log_R_a
    
    # Step 5: Calculate the final probability
    # P = C1 * p(z0, a1) + C2 * p(z0, a2)
    final_prob = C1 * p_z0_a1 + C2 * p_z0_a2

    print("Calculating the final probability P = C1*p(z0,a1) + C2*p(z0,a2):")
    print(f"p(z0, a1) = log({R}/{dist_z0_a1:.1f}) / log({R}) = {p_z0_a1:.4f}")
    print(f"p(z0, a2) = log({R}/{dist_z0_a2:.1f}) / log({R}) = {p_z0_a2:.4f}")
    print("\nFinal Equation:")
    print(f"P = {C1:.4f} * {p_z0_a1:.4f} + {C2:.4f} * {p_z0_a2:.4f}")
    print(f"P = {C:.4f} * ({p_z0_a1:.4f} + {p_z0_a2:.4f})")
    print(f"P = {final_prob:.4f}")
    print("-" * 30)

    # Round to three significant digits for the final answer
    print(f"The probability rounded to three significant digits is: {final_prob:.3g}")

solve_random_walk_problem()
<<<0.183>>>