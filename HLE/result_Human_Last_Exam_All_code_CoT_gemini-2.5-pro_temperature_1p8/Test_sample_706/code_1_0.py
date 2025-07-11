import math

def solve_random_walk_speed():
    """
    Calculates the asymptotic speed of a biased random walk on a random ladder graph.
    """
    # Probabilities from the problem description
    p_vertical_deleted = 1/2
    p_upper_horizontal_deleted = 1/3

    p_vertical_kept = 1 - p_vertical_deleted
    p_upper_horizontal_kept = 1 - p_upper_horizontal_deleted

    print("Step 1: Calculating the expected time to advance one unit horizontally.")
    print("----------------------------------------------------------------------")
    
    # E_0: Expected time starting from the lower rail (level 0)
    # The lower horizontal edge always exists, so it takes 1 step.
    E_0 = 1.0
    print(f"The expected time starting from the lower rail, E_0, is {E_0:.0f} step.")

    # E_1: Expected time starting from the upper rail (level 1)
    # This is a weighted average based on whether the upper horizontal edge exists.

    # First, calculate the expected time when the upper horizontal path is blocked.
    # The walker must backtrack k steps to find a vertical edge down.
    # Time taken T(k) = k_backward + 1_down + (k+1)_forward = 2k + 2 steps.
    # Prob of k backtracks P(k) = (p_vertical_deleted)^k * p_vertical_kept.
    # The expected time E_blocked = Sum[P(k) * T(k)] for k=0 to infinity.
    # This sum evaluates to 2 / p_vertical_kept.
    E_blocked = 2 / p_vertical_kept
    
    print(f"\nIf the upper path is blocked, the walker must find a vertical edge.")
    print(f"The probability a vertical edge exists is p_vertical_kept = {p_vertical_kept}.")
    print(f"The expected time to advance in this blocked case, E_blocked = 2 / p_vertical_kept = 2 / {p_vertical_kept} = {E_blocked:.1f} steps.")

    # E_1 = P(unblocked) * T(unblocked) + P(blocked) * T(blocked)
    E_1 = p_upper_horizontal_kept * 1 + p_upper_horizontal_deleted * E_blocked
    print(f"\nThe expected time from the upper rail, E_1, is calculated as:")
    print(f"E_1 = P(horizontal edge OK) * 1 + P(horizontal edge missing) * E_blocked")
    print(f"E_1 = {p_upper_horizontal_kept:.2f} * 1 + {p_upper_horizontal_deleted:.2f} * {E_blocked:.1f} = {E_1:.1f} steps.")
    print("\n")

    print("Step 2: Determining the walker's stationary vertical distribution.")
    print("----------------------------------------------------------------------")
    # In the c -> infinity limit, the walker never moves from the lower rail to the upper rail
    # because a path straight ahead is always available and infinitely preferred.
    # Thus, the lower rail is an absorbing state.
    # The stationary distribution (pi_0, pi_1) is (1, 0).
    pi_0 = 1.0
    pi_1 = 0.0
    print(f"The stationary distribution of being on the lower rail (pi_0) vs upper rail (pi_1) is:")
    print(f"(pi_0, pi_1) = ({pi_0:.1f}, {pi_1:.1f})")
    print("\n")

    print("Step 3: Calculating the final asymptotic speed v.")
    print("---------------------------------------------------")
    # Average time E[T] = pi_0 * E_0 + pi_1 * E_1
    E_T = pi_0 * E_0 + pi_1 * E_1
    print("The overall average time to advance one unit is E[T] = pi_0 * E_0 + pi_1 * E_1")
    print(f"E[T] = {pi_0:.1f} * {E_0:.1f} + {pi_1:.1f} * {E_1:.1f} = {E_T:.1f} steps.")

    # Asymptotic speed v = 1 / E[T]
    v = 1 / E_T
    print("\nThe asymptotic speed is v = 1 / E[T].")
    print(f"v = 1 / {E_T:.1f}")
    print(f"Final speed v = {v:.1f}")

solve_random_walk_speed()