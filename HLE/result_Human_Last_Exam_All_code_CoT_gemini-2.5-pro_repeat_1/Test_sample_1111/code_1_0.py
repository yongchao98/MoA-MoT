def solve_particle_problem():
    """
    This script explains the solution to the particle problem by analyzing
    the expected time T for different numbers of particles (k).
    """

    print("--- Analysis for k = 1 ---")
    print("If k=1, we have a single particle starting at a position x_1 > 0.")
    print("This particle performs a simple, unbiased random walk on the integers.")
    print("A standard result in probability theory is that for such a walk, the expected time to reach any specific state (like 0) is infinite.")
    print("Therefore, for k=1, E[T] is infinite.\n")

    print("--- Analysis for k = 2 ---")
    print("If k=2, we have particles at x_1 and x_2, with 0 < x_1 < x_2.")
    print("The first particle starts at x_1. Two things can happen:")
    print("1. The particle from x_1 reaches 0 before it reaches x_2.")
    print("2. The particle from x_1 reaches x_2 before it reaches 0.")
    print("\nThe expected time for one of these events to occur is finite. If event 1 occurs, the process stops, and the time taken is finite.")
    print("If event 2 occurs, the particle at x_2 is activated. At this moment, we have two active particles at position x_2.")
    print("Now, let's analyze the behavior of the minimum of these two active particles.")
    print("Let the two particles be at the same position 'x'. At each step, they move left or right independently.\n")
    
    print("We calculate the drift (expected change) of their minimum position:")
    # Probabilities of the two particles' movements (L for Left, R for Right)
    # P(L, L) = 1/2 * 1/2 = 1/4. New positions: (x-1, x-1). Minimum is x-1.
    # P(L, R) = 1/2 * 1/2 = 1/4. New positions: (x-1, x+1). Minimum is x-1.
    # P(R, L) = 1/2 * 1/2 = 1/4. New positions: (x+1, x-1). Minimum is x-1.
    # P(R, R) = 1/2 * 1/2 = 1/4. New positions: (x+1, x+1). Minimum is x+1.
    
    p_min_left = 0.75  # P(L,L) + P(L,R) + P(R,L) = 1/4 + 1/4 + 1/4
    p_min_right = 0.25 # P(R,R) = 1/4
    
    # Change in minimum position
    change_left = -1
    change_right = 1
    
    expected_change = (p_min_left * change_left) + (p_min_right * change_right)
    
    print("Let the two particles be at position x.")
    print(f"The probability that their minimum moves to x-1 is P(not both move right) = 1 - (1/2 * 1/2) = {p_min_left}")
    print(f"The probability that their minimum moves to x+1 is P(both move right) = 1/2 * 1/2 = {p_min_right}")
    print("\nThe expected change (drift) of the minimum position per step is:")
    print(f"E[change] = (P(min moves left) * (-1)) + (P(min moves right) * (+1))")
    print(f"E[change] = ({p_min_left} * {change_left}) + ({p_min_right} * {change_right}) = {expected_change}")
    
    print("\nSince the drift is negative, the minimum of the positions has a tendency to move towards 0.")
    print("A random walk with a negative drift has a finite expected time to reach 0 from any positive starting position.")
    print("Since all possible outcomes lead to a finite expected time, the total E[T] for k=2 is finite.")
    print("This conclusion does not depend on the specific initial positions x_1 and x_2.\n")

    print("--- Conclusion ---")
    print("E[T] is infinite for k=1 but finite for k=2.")
    minimal_k = 2
    print(f"The minimal value of k such that E[T] < infinity is {minimal_k}.")


solve_particle_problem()
<<<2>>>