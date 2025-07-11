import sympy
import itertools

def solve_random_walk_speed():
    """
    This script calculates the limit of the asymptotic speed v(c) of a biased
    random walk on a random ladder graph as the bias constant c tends to infinity.

    Method:
    1. The speed v(c) is expressed as a weighted average of speeds on the lower (0) and upper (1) levels:
       v(c) = pi_0(c) * E_0(c) + pi_1(c) * E_1(c)
       where pi_i is the stationary probability of being on level i, and E_i is the average speed on level i.

    2. The limits of E_0, E_1, pi_0, and pi_1 as c -> oo are calculated symbolically.

    3. The stationary probabilities are determined by balancing the transition rates between levels (R_01 and R_10).
       The ratio pi_1/pi_0 = R_01/R_10.

    4. All expectations (for speeds and rates) are computed by averaging over all possible local
       edge configurations, weighted by their probabilities.
    """
    
    # 1. Define symbolic variable c and model parameters
    c = sympy.Symbol('c', real=True, positive=True)
    exp_c = sympy.exp(c)
    exp_neg_c = sympy.exp(-c)

    p_upper_h = sympy.Rational(2, 3) # Probability of upper horizontal edge
    p_vertical = sympy.Rational(1, 2)  # Probability of vertical edge

    # Initialize symbolic expressions for expected speeds and transition rates
    E0, E1 = sympy.S(0), sympy.S(0)
    R01, R10 = sympy.S(0), sympy.S(0)

    # 2. Calculations for Level 0
    # The neighborhood of a vertex (n,0) depends on the existence of the vertical edge ((n,0),(n,1)).
    for v_n_exists in [0, 1]:
        p_config = p_vertical if v_n_exists else (1 - p_vertical)
        
        # Denominator for transition probabilities from (n,0)
        denom = exp_c + exp_neg_c + v_n_exists
        
        # Add contribution to the expected speed on level 0
        E0 += p_config * (exp_c - exp_neg_c) / denom
        
        # Add contribution to the transition rate from level 0 to 1
        R01 += p_config * v_n_exists / denom

    # 3. Calculations for Level 1
    # The neighborhood of (n,1) depends on three edges: U_n (right), U_nm1 (left), V_n (down).
    for u_n_exists, u_nm1_exists, v_n_exists in itertools.product([0, 1], repeat=3):
        p_config = (p_upper_h if u_n_exists else 1 - p_upper_h) * \
                   (p_upper_h if u_nm1_exists else 1 - p_upper_h) * \
                   (p_vertical if v_n_exists else 1 - p_vertical)
        
        denom = u_n_exists * exp_c + u_nm1_exists * exp_neg_c + v_n_exists
        # If denom is 0, the vertex is isolated, not part of the infinite component. Skip.
        if denom == 0:
            continue

        # Add contribution to the expected speed on level 1
        disp_num = u_n_exists * exp_c - u_nm1_exists * exp_neg_c
        E1 += p_config * disp_num / denom
        
        # Add contribution to the transition rate from level 1 to 0
        R10 += p_config * v_n_exists / denom

    # 4. Compute limits of all components as c -> infinity
    lim_E0 = sympy.limit(E0, c, sympy.oo)
    lim_E1 = sympy.limit(E1, c, sympy.oo)
    lim_R01 = sympy.limit(R01, c, sympy.oo)
    lim_R10 = sympy.limit(R10, c, sympy.oo)

    # 5. Determine the limit of the stationary probabilities
    # lim (pi1/pi0) = lim(R01)/lim(R10). We can swap limits due to continuity.
    lim_pi_ratio = lim_R01 / lim_R10
    lim_pi1 = lim_pi_ratio / (1 + lim_pi_ratio)
    lim_pi0 = 1 - lim_pi1
    
    # 6. Calculate the limit of the overall speed v(c)
    final_v = lim_pi0 * lim_E0 + lim_pi1 * lim_E1

    # 7. Print the results clearly, showing each component
    print("The final asymptotic speed v is calculated as: lim_{c->oo} v(c) = (lim pi_0) * (lim E_0) + (lim pi_1) * (lim E_1)")
    print("where pi_i are the stationary probabilities and E_i are the conditional average speeds.\n")
    print("The calculated limits of the components are:")
    print(f"lim E_0 = {lim_E0}")
    print(f"lim E_1 = {lim_E1}")
    print(f"lim R_01 (transition rate 0->1) = {lim_R01}")
    print(f"lim R_10 (transition rate 1->0) = {lim_R10}")
    print(f"lim pi_0 = {lim_pi0}")
    print(f"lim pi_1 = {lim_pi1}")
    print("\nSubstituting these values into the equation for the final speed:")
    print(f"{lim_pi0} * ({lim_E0}) + {lim_pi1} * ({lim_E1}) = {final_v}")

solve_random_walk_speed()