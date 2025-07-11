def solve_for_p():
    """
    This function calculates the price 'p' based on the provided economic model and data.
    """

    # Step 1: Define the known data and deduce basic parameters.
    # The joint probability distribution of (effort, output) pairs.
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # Identify distinct effort levels from the data.
    efforts = sorted(list(set(e for e, y in prob_dist.keys())))
    e_L, e_H = efforts[0], efforts[1]

    # Infer the states of the world 's' using the relationship y = s * e.
    # For e_H=22, y is 132 or 44, so s is 132/22=6 or 44/22=2.
    # This is consistent with e_L=10, where y is 60 or 20, giving s=6 or s=2.
    s_L, s_H = 2.0, 6.0

    print(f"From the data, we identified two effort levels: e_L = {e_L}, e_H = {e_H}")
    print(f"And two states of the world: s_L = {s_L}, s_H = {s_H}\n")


    # Step 2: Analyze employee's decision-making.
    # The employee's effort depends on the signal (theta_H or theta_L).
    # We associate high effort e_H with a "high" signal (theta_H) and low effort e_L with a "low" signal (theta_L).

    # Calculate the probability of observing each signal.
    # P(theta_H) = P(e=e_H)
    prob_theta_H = prob_dist[(e_H, s_H * e_H)] + prob_dist[(e_H, s_L * e_H)]
    # P(theta_L) = P(e=e_L)
    prob_theta_L = prob_dist[(e_L, s_H * e_L)] + prob_dist[(e_L, s_L * e_L)]

    # Calculate the employee's posterior beliefs P(s|theta).
    # P(s=s_H | theta=theta_H) = P(s=s_H, e=e_H) / P(e=e_H)
    prob_sH_given_thetaH = prob_dist[(e_H, s_H * e_H)] / prob_theta_H
    prob_sL_given_thetaH = prob_dist[(e_H, s_L * e_H)] / prob_theta_H

    # P(s=s_H | theta=theta_L) = P(s=s_H, e=e_L) / P(e=e_L)
    prob_sH_given_thetaL = prob_dist[(e_L, s_H * e_L)] / prob_theta_L
    prob_sL_given_thetaL = prob_dist[(e_L, s_L * e_L)] / prob_theta_L

    # Calculate the conditional expectation of the state given the signal.
    # E[s | theta_H]
    exp_s_given_thetaH = s_H * prob_sH_given_thetaH + s_L * prob_sL_given_thetaH
    # E[s | theta_L]
    exp_s_given_thetaL = s_H * prob_sH_given_thetaL + s_L * prob_sL_given_thetaL

    print(f"Calculated conditional expectations of the state:")
    print(f"E[s | signal=high] = {exp_s_given_thetaH}")
    print(f"E[s | signal=low] = {exp_s_given_thetaL}\n")

    # Step 3: Determine the contract parameter beta.
    # The employee's optimal effort is e = beta * E[s|theta].
    # We can use the high signal case to find beta.
    beta = e_H / exp_s_given_thetaH
    # We can check this with the low signal case: beta should be e_L / exp_s_given_thetaL
    
    print(f"Using the employee's optimal effort rule e = beta * E[s|theta], we find beta:")
    print(f"{e_H} = beta * {exp_s_given_thetaH}  =>  beta = {beta}\n")

    # Step 4: Solve the firm's problem.
    # The firm chooses beta to maximize profit. The first-order condition for this problem
    # shows that the optimal beta is equal to the price, p.
    # This is because the firm balances the marginal revenue from higher effort (induced by a higher beta)
    # against the marginal cost of higher effort (which the firm must compensate the employee for).
    # This leads to the condition p = beta.

    p = beta

    print("The firm's profit maximization leads to the condition p = beta.")
    print("Therefore, the value of p is:")
    print(f"p = {p}")

solve_for_p()
>>>4