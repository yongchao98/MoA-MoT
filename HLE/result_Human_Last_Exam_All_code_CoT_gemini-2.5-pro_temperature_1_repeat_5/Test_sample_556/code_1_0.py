import collections

def solve_for_p():
    """
    This function solves for the price 'p' based on the provided economic model and data.
    """
    # Step 1: Deconstruct the Data
    # Given probability distribution P(e, y)
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # From the data, we can identify the high and low effort levels.
    e_H = 22
    e_L = 10

    # From y = s * e, we can identify the states of the world s_H and s_L.
    # Using the first data point with high effort: y = 132, e = 22 => s = 132 / 22 = 6
    s_H = 132 / e_H
    # Using the second data point with high effort: y = 44, e = 22 => s = 44 / 22 = 2
    s_L = 44 / e_H

    print(f"Step 1: Identified model parameters from data")
    print(f"High effort (e_H): {e_H}")
    print(f"Low effort (e_L): {e_L}")
    print(f"High state (s_H): {s_H}")
    print(f"Low state (s_L): {s_L}\n")

    # Step 2: Determine Incentive Parameter β
    # Map the given probabilities to the joint probabilities of signal (θ) and state (s)
    # P(e=e_H, y=s_H*e_H) = P(θ=s_H, s=s_H)
    # P(e=e_H, y=s_L*e_H) = P(θ=s_H, s=s_L)
    # P(e=e_L, y=s_H*e_L) = P(θ=s_L, s=s_H)
    # P(e=e_L, y=s_L*e_L) = P(θ=s_L, s=s_L)
    P_thetaH_sH = prob_dist[(e_H, s_H * e_H)]
    P_thetaH_sL = prob_dist[(e_H, s_L * e_H)]
    P_thetaL_sH = prob_dist[(e_L, s_H * e_L)]
    P_thetaL_sL = prob_dist[(e_L, s_L * e_L)]

    # Calculate marginal probabilities of signals
    P_thetaH = P_thetaH_sH + P_thetaH_sL
    P_thetaL = P_thetaL_sH + P_thetaL_sL

    # Calculate conditional probabilities P(s|θ) using Bayes' rule
    P_sH_given_thetaH = P_thetaH_sH / P_thetaH
    P_sL_given_thetaH = P_thetaH_sL / P_thetaH
    P_sH_given_thetaL = P_thetaL_sH / P_thetaL
    P_sL_given_thetaL = P_thetaL_sL / P_thetaL

    # Calculate expected state E[s|θ] for each signal
    E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
    E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL

    # From the employee's effort choice e = β * E[s|θ], we can find β
    beta = e_H / E_s_given_thetaH
    
    print(f"Step 2: Determined incentive parameter β")
    print(f"Expected state given high signal E[s|θ=s_H]: {E_s_given_thetaH}")
    print(f"Expected state given low signal E[s|θ=s_L]: {E_s_given_thetaL}")
    print(f"From e_H = β * E[s|θ=s_H] => {e_H} = β * {E_s_given_thetaH}")
    print(f"Solved β = {beta:.2f}\n")

    # Step 3: Determine Fixed Payment α
    # Employee's expected utility, E[U], is set to 0 by the participation constraint.
    # E[U|θ] = α + (β^2 * E[s|θ]^2) / 2
    E_U_given_thetaH = (beta**2 * E_s_given_thetaH**2) / 2
    E_U_given_thetaL = (beta**2 * E_s_given_thetaL**2) / 2
    # E[U] = α + P(θ=s_H)*E[U|θ=s_H] + P(θ=s_L)*E[U|θ=s_L]
    E_U_without_alpha = P_thetaH * E_U_given_thetaH + P_thetaL * E_U_given_thetaL
    # E[U] = α + E_U_without_alpha = 0 => α = -E_U_without_alpha
    alpha = -E_U_without_alpha

    print(f"Step 3: Determined fixed payment α")
    print(f"Employee's expected utility is E[U] = α + {E_U_without_alpha:.2f}")
    print(f"The participation constraint implies E[U] = 0, so α = {alpha:.2f}\n")

    # Step 4: Determine Price p
    # Calculate expected output E[y]
    E_y = sum(y * prob for (e, y), prob in prob_dist.items())

    # The firm's expected profit is E[Π] = (p - β)E[y] - α.
    # Assuming a competitive market, E[Π] = 0.
    # (p - β)E[y] - α = 0 => p - β = α / E[y] => p = β + α / E[y]
    p = beta + alpha / E_y

    print("Step 4: Determine price p using the zero-profit condition")
    print(f"Expected output E[y]: {E_y:.2f}")
    print("The firm's expected profit is E[Π] = (p - β) * E[y] - α.")
    print("Setting E[Π] = 0 gives the equation: (p - β) * E[y] = α")
    print(f"Substituting the calculated values:")
    print(f"(p - {beta:.2f}) * {E_y:.2f} = {alpha:.2f}")
    
    val = alpha / E_y
    print(f"p - {beta:.2f} = {alpha:.2f} / {E_y:.2f} = {val:.2f}")
    print(f"p = {beta:.2f} + ({val:.2f})")
    print(f"p = {p:.2f}\n")

    # Final Answer
    print(f"The value of p is {p}.")
    return p

if __name__ == '__main__':
    final_p = solve_for_p()
    # The final answer is requested in a specific format.
    # print(f"<<<{final_p}>>>") # This is for the platform, not for user copy-paste.
