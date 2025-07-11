import collections

def solve_for_p():
    """
    This function solves for the price p based on the provided economic model and data.
    """
    # Step 1: Deconstruct the Data
    # Given probability distribution P(e, y)
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375,
    }

    # Extract unique effort levels and identify high/low effort
    effort_levels = sorted(list(set(e for e, y in prob_dist.keys())), reverse=True)
    e_H, e_L = effort_levels[0], effort_levels[1]

    # Determine states of the world s = y / e
    states = sorted(list(set(y / e for e, y in prob_dist.keys())), reverse=True)
    s_H, s_L = states[0], states[1]

    print("Step 1: Extracting model parameters from data.")
    print(f"High effort level (e_H): {e_H}")
    print(f"Low effort level (e_L): {e_L}")
    print(f"High state of the world (s_H): {s_H}")
    print(f"Low state of the world (s_L): {s_L}")
    print("-" * 30)

    # Step 2: Determine model probabilities
    # We associate high effort e_H with signal theta_H and low effort e_L with signal theta_L.
    # The probability P(e,y) corresponds to P(theta, s)
    # P(theta=s_H, s=s_H)
    p_thetaH_sH = prob_dist[(e_H, e_H * s_H)]
    # P(theta=s_H, s=s_L)
    p_thetaH_sL = prob_dist[(e_H, e_H * s_L)]
    # P(theta=s_L, s=s_H)
    p_thetaL_sH = prob_dist[(e_L, e_L * s_H)]
    # P(theta=s_L, s=s_L)
    p_thetaL_sL = prob_dist[(e_L, e_L * s_L)]

    # Calculate marginal probability of receiving signal theta_H
    p_thetaH = p_thetaH_sH + p_thetaH_sL

    # Calculate the conditional expectation of the state given the high signal, E[s | theta_H]
    # E[s|theta_H] = s_H * P(s_H|theta_H) + s_L * P(s_L|theta_H)
    #              = s_H * (P(s_H, theta_H)/P(theta_H)) + s_L * (P(s_L, theta_H)/P(theta_H))
    E_s_given_thetaH = (s_H * p_thetaH_sH + s_L * p_thetaH_sL) / p_thetaH

    print("Step 2: Calculating conditional expectation of the state.")
    print(f"Probability of receiving high signal P(theta_H): {p_thetaH}")
    print(f"Expected state given high signal E[s|theta_H]: {E_s_given_thetaH}")
    print("-" * 30)
    
    # Step 3: Find the incentive parameter beta
    # The employee's optimal effort is e(theta) = beta * E[s|theta].
    # Using the high effort and high signal case: e_H = beta * E[s|theta_H]
    beta = e_H / E_s_given_thetaH

    print("Step 3: Calculating the incentive parameter beta.")
    print(f"From the employee's effort choice rule, e_H = beta * E[s|theta_H].")
    print(f"So, beta = e_H / E[s|theta_H] = {e_H} / {E_s_given_thetaH} = {beta}")
    print("-" * 30)

    # Step 4 & 5: Analyze the Firm's Problem and Solve for p
    # The firm chooses beta to maximize its profit. The first-order condition for this maximization
    # shows that the optimal choice for the incentive parameter is beta = p.
    # The observed data results from this optimal choice.
    p = beta
    
    print("Step 4 & 5: Determining the price p.")
    print("The firm chooses beta to maximize its profit.")
    print("The profit-maximizing condition is for the firm to set the incentive parameter beta equal to the price p.")
    print("Since the observed data comes from this optimal choice, we have p = beta.")
    print("\nFinal Calculation:")
    print(f"p = {p}")

solve_for_p()