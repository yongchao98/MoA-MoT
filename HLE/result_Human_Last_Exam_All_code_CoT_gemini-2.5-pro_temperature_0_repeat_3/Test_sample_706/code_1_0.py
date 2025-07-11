import math

def solve_random_walk():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.
    """
    # Step 1: Define the probabilities of edge existence from the problem statement.
    # Probability a vertical edge exists
    p_v = 1 - 1/2
    # Probability an upper horizontal edge exists
    p_h = 1 - 1/3
    
    print("Step 1: Analyzing the stationary distribution pi_i(c) for rail i as c -> infinity.")
    
    # Step 2: Analyze the transition rates between rails for large c.
    # Rate from bottom to top (0 -> 1) is proportional to exp(-c) because moving up (weight 1)
    # competes with moving right (weight exp(c)).
    # k_01 is proportional to (1/2) * 1 / (exp(c) + exp(-c) + 1) ~ exp(-c)
    print("The transition rate from bottom to top, k_01(c), is proportional to exp(-c).")
    
    # Rate from top to bottom (1 -> 0) happens when a horizontal path is blocked
    # and a vertical path exists.
    # This occurs with a probability independent of c for large c.
    # P(horizontal blocked) = 1 - p_h = 1/3
    # P(vertical exists) = p_v = 1/2
    # When this happens, the move down (weight 1) is overwhelmingly preferred to moving left (weight exp(-c)).
    k10_prob = (1 - p_h) * p_v
    print(f"The transition rate from top to bottom, k_10(c), approaches a constant for large c.")
    print(f"This rate is proportional to P(horizontal blocked) * P(vertical exists) = {1-p_h:.2f} * {p_v:.2f} = {k10_prob:.2f}.")

    # In equilibrium, pi_0 * k_01 = pi_1 * k_10.
    # So, pi_1 / pi_0 = k_01 / k_10, which is proportional to exp(-c).
    # As c -> infinity, this ratio goes to 0.
    print("\nThe ratio pi_1(c)/pi_0(c) is proportional to exp(-c), which tends to 0 as c -> infinity.")
    
    # Since pi_0 + pi_1 = 1, if their ratio is 0, then pi_0 must be 1 and pi_1 must be 0.
    lim_pi_0 = 1
    lim_pi_1 = 0
    print(f"Therefore, in the limit c -> infinity, the stationary probabilities are:")
    print(f"lim pi_0(c) = {lim_pi_0}")
    print(f"lim pi_1(c) = {lim_pi_1}")
    print("This means the walker is almost certainly on the bottom rail for very large c.")
    
    print("\nStep 2: Analyzing the speed v_i(c) on each rail as c -> infinity.")
    
    # Step 3: Analyze the speed on the bottom rail.
    # The bottom rail is always intact. The walker almost always moves right.
    lim_v_0 = 1
    print(f"On the bottom rail, rightward edges always exist. The speed v_0(c) approaches {lim_v_0} as c -> infinity.")
    
    # Step 4: Analyze the speed on the top rail.
    # The top rail has traps. A trap occurs if a rightward horizontal edge is missing,
    # but the preceding one existed, and the vertical edge is also missing.
    # P(trap) = P(H_exist) * P(H_missing) * P(V_missing)
    prob_trap_config = p_h * (1 - p_h) * (1 - p_v)
    # The time spent in a trap is proportional to exp(c). This makes the speed on the top rail go to 0.
    lim_v_1 = 0
    print(f"On the top rail, traps exist with probability {p_h:.2f} * {1-p_h:.2f} * {1-p_v:.2f} = {prob_trap_config:.3f} per site.")
    print(f"The time spent in these traps is exponential in c, so the speed v_1(c) approaches {lim_v_1} as c -> infinity.")
    
    # Step 5: Combine the results.
    # v(c) = pi_0(c)*v_0(c) + pi_1(c)*v_1(c)
    # lim v(c) = (lim pi_0)*(lim v_0) + (lim pi_1)*(lim v_1)
    final_speed = lim_pi_0 * lim_v_0 + lim_pi_1 * lim_v_1
    
    print("\nStep 3: Combining the results to find the final speed limit.")
    print(f"The final speed is lim v(c) = (lim pi_0(c)) * (lim v_0(c)) + (lim pi_1(c)) * (lim v_1(c)).")
    print(f"The final equation is: {lim_pi_0} * {lim_v_0} + {lim_pi_1} * {lim_v_1}")
    print(f"The result is {final_speed}.")
    
    return final_speed

solve_random_walk()
<<<1>>>