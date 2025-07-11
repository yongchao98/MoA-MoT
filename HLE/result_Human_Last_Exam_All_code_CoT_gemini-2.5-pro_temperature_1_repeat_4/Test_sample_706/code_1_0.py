import math

def solve():
    """
    Calculates the limiting speed v(c) as c -> infinity.
    The logic follows the explained steps.
    """

    # Probability that an upper horizontal edge exists
    p_u = 1 - 1/3

    # Probability that a vertical edge exists
    p_v = 1 - 1/2

    # Step 1: Speed on the lower rail (v0)
    # As c -> infinity, the walker on the lower rail always moves right.
    # Displacement is 1 per step. Time is 1 per step.
    v0 = 1.0
    print(f"The asymptotic speed on the lower rail (v0) is {v0}.")

    # Step 2: Speed on the upper rail (v1)
    # This is calculated based on the 'intelligent walker' model for c -> infinity.
    
    # Expected length of a forward run on the upper rail
    E_N = p_u / (1 - p_u)
    
    # Expected length of a backward search for a downward exit
    E_L = 1 / p_v
    
    # Expected displacement for a full cycle on the upper rail
    # A cycle consists of a forward run followed by either an escape or a backtrack-and-escape.
    # The displacement is E[N] from the run, and with probability (1-p_v) an additional -E[L] from backtracking.
    expected_displacement_v1 = E_N - (1 - p_v) * E_L
    
    # Expected time for a full cycle on the upper rail
    # Time is E[N] for the run.
    # With prob p_v, 1 extra step to go down.
    # With prob (1-p_v), E[L] steps to backtrack + 1 step to go down.
    expected_time_v1 = E_N + p_v * 1 + (1 - p_v) * (E_L + 1)
    
    # The speed on the upper rail is the ratio of expected displacement to expected time.
    v1 = expected_displacement_v1 / expected_time_v1
    print(f"The asymptotic speed on the upper rail (v1) is {expected_displacement_v1} / {expected_time_v1} = {v1}.")

    # Step 3: Combine speeds for the final answer
    # As c -> infinity, the probability of moving from the lower rail to the upper rail
    # approaches 0 (proportional to exp(-c)).
    # Therefore, the stationary probability of being on the upper rail, pi_1, goes to 0.
    # The walker spends effectively 100% of its time on the lower rail.
    pi_0 = 1.0
    pi_1 = 0.0
    
    final_speed = pi_0 * v0 + pi_1 * v1
    
    print(f"\nThe final limiting speed is (pi_0 * v0) + (pi_1 * v1)")
    print(f"= ({pi_0} * {v0}) + ({pi_1} * {v1})")
    print(f"= {final_speed}")

solve()