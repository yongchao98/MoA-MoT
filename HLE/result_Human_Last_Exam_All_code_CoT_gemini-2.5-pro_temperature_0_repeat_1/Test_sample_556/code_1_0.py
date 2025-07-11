import collections

def solve_price():
    """
    Calculates the price p based on the provided probability distribution.
    """
    # Step 1: Define the given joint probability distribution P(e, y)
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # Step 2: Determine the states of the world 's' and effort levels 'e' from the data.
    # The relation is y = s * e.
    # We can find the unique values for s and e from the keys of the dictionary.
    efforts = sorted(list(set(e for e, y in prob_dist.keys())), reverse=True)
    e_H, e_L = efforts[0], efforts[1] # High and low effort

    # Calculate s for a given effort level
    s_values_for_e_H = [y / e_H for e, y in prob_dist.keys() if e == e_H]
    states = sorted(list(set(s_values_for_e_H)), reverse=True)
    s_H, s_L = states[0], states[1] # High and low state

    print(f"From the data, we deduce the following:")
    print(f"High effort e_H = {e_H}, Low effort e_L = {e_L}")
    print(f"High state s_H = {s_H}, Low state s_L = {s_L}\n")

    # Step 3: Map the (e, y) probabilities to (theta, s) probabilities.
    # High effort (e_H) corresponds to the high signal (theta_H).
    # Low effort (e_L) corresponds to the low signal (theta_L).
    p_thetaH_sH = prob_dist[(e_H, s_H * e_H)]
    p_thetaH_sL = prob_dist[(e_H, s_L * e_H)]
    p_thetaL_sH = prob_dist[(e_L, s_H * e_L)]
    p_thetaL_sL = prob_dist[(e_L, s_L * e_L)]

    # Step 4: Calculate the marginal probability of receiving the high signal, P(theta_H).
    p_thetaH = p_thetaH_sH + p_thetaH_sL
    print(f"The probability of the agent receiving the high signal P(θ=s_H) is: {p_thetaH}\n")

    # Step 5: Calculate the conditional expectation E[s | theta_H].
    # E[s|theta_H] = s_H * P(s_H|theta_H) + s_L * P(s_L|theta_H)
    # P(s|theta_H) = P(s, theta_H) / P(theta_H)
    p_sH_given_thetaH = p_thetaH_sH / p_thetaH
    p_sL_given_thetaH = p_thetaH_sL / p_thetaH
    E_s_given_thetaH = s_H * p_sH_given_thetaH + s_L * p_sL_given_thetaH
    
    print(f"The agent's conditional expectation of the state, given the high signal, is:")
    print(f"E[s|θ=s_H] = {s_H} * P(s=s_H|θ=s_H) + {s_L} * P(s=s_L|θ=s_H)")
    print(f"E[s|θ=s_H] = {s_H} * {p_sH_given_thetaH:.4f} + {s_L} * {p_sL_given_thetaH:.4f} = {E_s_given_thetaH}\n")

    # Step 6: Solve for beta using the agent's optimal effort rule: e_H = beta * E[s | theta_H].
    # The firm's profit maximization implies p = beta.
    beta = e_H / E_s_given_thetaH
    p = beta

    print("The agent's optimal effort rule is e = β * E[s|θ].")
    print("The firm's profit maximization implies that the price p is equal to β.")
    print("We can solve for β using the high effort and high signal case:")
    print(f"β = e_H / E[s|θ=s_H]")
    print(f"β = {e_H} / {E_s_given_thetaH}")
    print(f"p = β = {beta}\n")
    
    print("The final value of p is:")
    print(int(p))

solve_price()