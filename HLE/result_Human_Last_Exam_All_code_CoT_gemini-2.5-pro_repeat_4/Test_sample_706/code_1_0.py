from fractions import Fraction

def solve():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.
    """

    # Probabilities of edges existing
    p_h = Fraction(2, 3)  # Probability of upper horizontal edge
    p_v = Fraction(1, 2)  # Probability of vertical edge

    # Step 1: Calculate v_0_lim, the speed on the lower level as c -> inf.
    # On the lower level, the horizontal edge to the right always exists.
    # The weight for moving right is e^c, for moving left e^{-c}, and for moving up is 1.
    # As c -> inf, the probability of moving right approaches 1.
    # The displacement is +1 in one time step.
    v_0_lim = 1
    print(f"Limit of speed on lower level (v_0_lim): {v_0_lim}")

    # Step 2: Calculate v_1_lim, the speed on the upper level as c -> inf.
    # This is the expected displacement, averaged over all local edge configurations.
    # Let H_n be the right edge, H'_{n-1} the left edge, V_n the vertical edge.
    
    # Case 1: H_n exists (prob p_h). Limiting displacement is +1.
    v1_case1 = 1 * p_h
    
    # Case 2: H_n is missing, V_n exists (prob (1-p_h)*p_v). Limiting displacement is 0.
    v1_case2 = 0 * (1 - p_h) * p_v
    
    # Case 3: H_n missing, V_n missing, H'_{n-1} exists (prob (1-p_h)*(1-p_v)*p_h). Limiting displacement is -1.
    v1_case3 = -1 * (1 - p_h) * (1 - p_v) * p_h
    
    # A vertex can be isolated from its neighbours if H_n, H'_{n-1}, V_n are all missing.
    # The walk is on the infinite component, so we ignore these isolated vertices.
    prob_isolated = (1 - p_h) * (1 - p_h) * (1 - p_v)
    
    # The expected displacement is the sum of contributions, normalized by the probability of being connected.
    v_1_lim_unnormalized = v1_case1 + v1_case2 + v1_case3
    normalization_factor = 1 / (1 - prob_isolated)
    
    v_1_lim = v_1_lim_unnormalized * normalization_factor
    print(f"Limit of speed on upper level (v_1_lim): {v_1_lim.numerator}/{v_1_lim.denominator}")

    # Step 3: Analyze the stationary distribution pi = (pi_0, pi_1).
    # From detailed balance, pi_1/pi_0 = P(0->1)/P(1->0). We need the limits of these probabilities.
    # P(0->1|c) is proportional to e^{-c}, so its limit is 0.
    # Let's find the limit of P(1->0|c). It's non-zero only when H_n is missing and V_n exists.
    # This occurs with probability (1 - p_h) * p_v. In this case, the jump probability -> 1.
    p_1_to_0_unnormalized = 1 * (1 - p_h) * p_v
    p_1_to_0_lim = p_1_to_0_unnormalized * normalization_factor
    
    print(f"Limit of transition probability from level 1 to 0 (P_1_to_0_lim): {p_1_to_0_lim.numerator}/{p_1_to_0_lim.denominator}")

    # Since P(0->1) -> 0 and P(1->0) -> constant > 0, we have lim pi_1 = 0.
    pi_1_lim = 0
    pi_0_lim = 1
    print(f"Limit of stationary probability on upper level (pi_1_lim): {pi_1_lim}")
    print(f"Limit of stationary probability on lower level (pi_0_lim): {pi_0_lim}")

    # Step 4: Combine the results to get the final speed limit.
    # v_inf = pi_0_lim * v_0_lim + pi_1_lim * v_1_lim
    final_v_lim = pi_0_lim * v_0_lim + pi_1_lim * v_1_lim
    
    print(f"\nThe final equation is: v_inf = pi_0_lim * v_0_lim + pi_1_lim * v_1_lim")
    print(f"v_inf = {pi_0_lim} * {v_0_lim} + {pi_1_lim} * {v_1_lim.numerator}/{v_1_lim.denominator}")
    print(f"Final Result (v_inf): {final_v_lim}")

solve()