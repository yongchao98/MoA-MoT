import numpy as np

def calculate_speed_limit():
    """
    Calculates the asymptotic speed v(c) by analyzing the stationary state
    of the system for a very large constant c.
    """
    
    # We choose a large value for c to approximate the limit c -> infinity
    c = 100.0
    
    # Probabilities of edges existing
    p_v = 1 / 2  # Vertical edge
    p_u = 2 / 3  # Upper horizontal edge

    ec = np.exp(c)
    emc = np.exp(-c)

    # In the stationary state, probability flows between levels must balance.
    # Flow_up = rho_0 * T_01, Flow_down = rho_1 * T_10
    
    # T_01: Average transition rate from level 0 to 1.
    # Averages over the two cases for the vertical edge V_n at a site n.
    prob_up_if_v_exists = 1 / (ec + emc + 1)
    prob_up_if_v_absent = 0
    T_01 = p_v * prob_up_if_v_exists + (1 - p_v) * prob_up_if_v_absent

    # T_10: Average transition rate from level 1 to 0.
    # Averages over all 8 local graph configurations for (V_n, H_n, H_{n-1}).
    prob_down = 0
    configs = [(v, h_n, h_nm1) for v in [0, 1] for h_n in [0, 1] for h_nm1 in [0, 1]]
    for v, h_n, h_nm1 in configs:
        p_config = (p_v if v == 1 else 1 - p_v) * \
                   (p_u if h_n == 1 else 1 - p_u) * \
                   (p_u if h_nm1 == 1 else 1 - p_u)
        
        if v == 1: # No downward move if vertical edge is absent
            denominator = h_n * ec + h_nm1 * emc + v
            prob_down += p_config * (v / denominator)
    T_10 = prob_down

    # Solve rho_0 * T_01 = rho_1 * T_10 and rho_0 + rho_1 = 1
    # This gives rho_1 = T_01 / (T_10 + T_01)
    if T_10 + T_01 == 0:
      rho_1 = 0.0
    else:
      rho_1 = T_01 / (T_10 + T_01)
    rho_0 = 1 - rho_1

    # Calculate v0, the average speed on level 0.
    # v0 = E[P(right|0) - P(left|0)]
    v0 = (p_v * (ec-emc)/(ec+emc+1)) + ((1-p_v) * (ec-emc)/(ec+emc))

    # Calculate v1, the average speed on level 1.
    # v1 = E[P(right|1) - P(left|1)]
    v1 = 0
    for v, h_n, h_nm1 in configs:
        p_config = (p_v if v == 1 else 1 - p_v) * \
                   (p_u if h_n == 1 else 1 - p_u) * \
                   (p_u if h_nm1 == 1 else 1 - p_u)
        denominator = h_n * ec + h_nm1 * emc + v
        if denominator > 0:
            v1 += p_config * (h_n * ec - h_nm1 * emc) / denominator

    # The total asymptotic speed v(c) = rho_0 * v0 + rho_1 * v1
    v = rho_0 * v0 + rho_1 * v1

    print(f"For a large constant c = {c}:")
    print(f"The stationary probability of being on the lower rail is rho_0 = {rho_0}")
    print(f"The stationary probability of being on the upper rail is rho_1 = {rho_1}")
    print(f"The average speed on the lower rail is v0 = {v0}")
    print(f"The average speed on the upper rail is v1 = {v1}")
    print("\nThe total asymptotic speed v(c) is calculated as rho_0 * v0 + rho_1 * v1:")
    print(f"v({c}) = {rho_0:.5f} * {v0:.5f} + {rho_1:.5e} * {v1:.5f}")
    print(f"v({c}) = {v}")
    
    # The final answer is the limit of v(c) as c -> infinity
    final_answer = round(v)
    print(f"\nAs c approaches infinity, v(c) converges to {final_answer}.")
    return final_answer
    
calculate_speed_limit()