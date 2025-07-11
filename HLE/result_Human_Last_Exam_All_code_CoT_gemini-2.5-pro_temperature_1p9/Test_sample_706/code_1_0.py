import sympy

def solve_speed_limit():
    """
    Calculates the limit of the asymptotic speed v(c) as c -> infinity.
    The method follows the analytical derivation explained above.
    """

    # --- Step 1: Define parameters ---
    # Probability that an upper horizontal edge exists
    p_u = sympy.Rational(2, 3)
    # Probability that a vertical edge exists
    p_v = sympy.Rational(1, 2)
    
    print("--- Analytical Calculation of lim_{c->inf} v(c) ---")
    print(f"Probability of upper horizontal edge: {p_u}")
    print(f"Probability of vertical edge: {p_v}\n")

    # --- Step 2: Calculate lim v_0(c) ---
    # As c -> inf, displacement from the lower rail always tends to 1.
    v0_lim = 1
    print(f"Step 2: Speed on the lower rail, v_0(inf)")
    print(f"The limit of the expected displacement from the lower rail is 1.")
    print(f"lim v_0(c) = {v0_lim}\n")

    # --- Step 3: Calculate lim v_1(c) ---
    # As c -> inf, displacement from upper rail depends on neighbors.
    # Case 1: Forward upper edge exists. Displacement limit = 1.
    prob_v1_forward = p_u
    disp_v1_forward = 1

    # Case 2: Forward upper edge is missing, have to go back. Disp limit = -1.
    # This happens if U_{n+1}=0, U_n=1, V_n=0.
    prob_v1_backward = (1 - p_u) * p_u * (1 - p_v)
    disp_v1_backward = -1

    # Case 3: Otherwise, displacement limit is 0.
    v1_lim = prob_v1_forward * disp_v1_forward + prob_v1_backward * disp_v1_backward
    
    print("Step 3: Speed on the upper rail, v_1(inf)")
    print(f"This is E[lim T], where T is the displacement from the upper rail.")
    print(f"Term 1: Move Forward => P = {prob_v1_forward}, Disp = {disp_v1_forward}")
    print(f"Term 2: Move Backward => P = {prob_v1_backward}, Disp = {disp_v1_backward}")
    print(f"Term 3: Move Vertically or Trapped => P = 1 - ({prob_v1_forward} + {prob_v1_backward}), Disp = 0")
    print(f"lim v_1(c) = ({prob_v1_forward})*({disp_v1_forward}) + ({prob_v1_backward})*({disp_v1_backward}) + ... = {v1_lim}\n")

    # --- Step 4: Calculate lim pi_0 and pi_1 ---
    # From the flux balance equation pi_0 * E[P_01] = pi_1 * E[P_10],
    # E[P_01] is proportional to exp(-c) while E[P_10] tends to a constant.
    # Therefore, pi_1 / pi_0 tends to 0.
    pi0_lim = 1
    pi1_lim = 0
    print("Step 4: Stationary distribution on rails, pi_0(inf), pi_1(inf)")
    print("The flux balance shows that the probability of being on the upper rail vanishes.")
    print(f"lim pi_0(c) = {pi0_lim}")
    print(f"lim pi_1(c) = {pi1_lim}\n")
    
    # --- Step 5: Final Result ---
    v_lim = pi0_lim * v0_lim + pi1_lim * v1_lim
    print("Step 5: Final Speed Limit Calculation")
    print("v(inf) = pi_0(inf) * v_0(inf) + pi_1(inf) * v_1(inf)")
    print(f"v(inf) = ({pi0_lim}) * ({v0_lim}) + ({pi1_lim}) * ({v1_lim})")
    print(f"v(inf) = {v_lim}")
    
    return v_lim

if __name__ == '__main__':
    final_speed = solve_speed_limit()
    print(f"\nFinal Answer: The limit of the asymptotic speed is {final_speed}.")
