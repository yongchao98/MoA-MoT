import collections

def solve_price():
    """
    This function solves for the price 'p' based on the provided economic model and data.
    """
    # Step 1: Deconstruct the Data
    # Given joint probability distribution P(e, y)
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # Identify effort levels. We assume higher effort corresponds to the high signal.
    e_H = 22
    e_L = 10
    print(f"Effort levels are e_H = {e_H} and e_L = {e_L}.")

    # From y = s * e, we can determine the states of the world, s.
    # For e = 22, y can be 132 or 44. So s can be 132/22=6 or 44/22=2.
    # For e = 10, y can be 60 or 20. So s can be 60/10=6 or 20/10=2.
    s_H = 6
    s_L = 2
    print(f"States of the world are s_H = {s_H} and s_L = {s_L}.")

    # Calculate probabilities of receiving each signal theta.
    # P(theta_H) = P(e=e_H)
    p_theta_H = prob_dist[(22, 132)] + prob_dist[(22, 44)]
    # P(theta_L) = P(e=e_L)
    p_theta_L = prob_dist[(10, 60)] + prob_dist[(10, 20)]
    print(f"Probability of high signal P(theta_H) = {p_theta_H:.4f}")
    print(f"Probability of low signal P(theta_L) = {p_theta_L:.4f}")

    # Step 2: Employee Behavior (already described in the plan)
    print("\nEmployee chooses effort 'e' to maximize E[w - e^2/2 | theta].")
    print("This leads to the optimality condition: e = beta * E[s | theta].")

    # Step 3: Calculate Conditional Expectations E[s|theta]
    # P(s, theta) corresponds to the given P(e, y)
    # P(s_H, theta_H) corresponds to P(e_H, y=s_H*e_H) = P(22, 132)
    p_sH_thetaH = prob_dist[(22, 132)]
    # P(s_L, theta_H) corresponds to P(e_H, y=s_L*e_H) = P(22, 44)
    p_sL_thetaH = prob_dist[(22, 44)]

    # P(s_H, theta_L) corresponds to P(e_L, y=s_H*e_L) = P(10, 60)
    p_sH_thetaL = prob_dist[(10, 60)]
    # P(s_L, theta_L) corresponds to P(e_L, y=s_L*e_L) = P(10, 20)
    p_sL_thetaL = prob_dist[(10, 20)]

    # E[s | theta_H] = s_H * P(s_H | theta_H) + s_L * P(s_L | theta_H)
    # where P(s | theta) = P(s, theta) / P(theta)
    E_s_given_theta_H = s_H * (p_sH_thetaH / p_theta_H) + s_L * (p_sL_thetaH / p_theta_H)
    
    # E[s | theta_L] = s_H * P(s_H | theta_L) + s_L * P(s_L | theta_L)
    E_s_given_theta_L = s_H * (p_sH_thetaL / p_theta_L) + s_L * (p_sL_thetaL / p_theta_L)

    print(f"\nCalculated conditional expectation E[s | theta_H] = {E_s_given_theta_H}")
    print(f"Calculated conditional expectation E[s | theta_L] = {E_s_given_theta_L}")

    # Step 4: Solve for beta
    # Using the equation for the high signal: e_H = beta * E[s | theta_H]
    beta = e_H / E_s_given_theta_H
    print(f"\nUsing the optimality condition for the high signal:")
    print(f"{e_H} = beta * {E_s_given_theta_H}")
    print(f"Solving for beta, we get beta = {int(beta)}.")

    # We can check this with the low signal equation: e_L = beta * E[s | theta_L]
    beta_check = e_L / E_s_given_theta_L
    print(f"Checking with the low signal: {e_L} = {int(beta_check)} * {E_s_given_theta_L}, which is consistent.")

    # Step 5 & 6: Firm Behavior and Solving for p
    # For a linear contract w = alpha + beta*y and quadratic effort cost e^2/2,
    # the firm's profit maximization problem with respect to beta leads to the
    # first-order condition: p = beta.
    p = beta
    print("\nThe firm chooses beta to maximize expected profit.")
    print("The first-order condition for this model is p = beta.")
    print("\nTherefore, the final equation is:")
    print(f"p = {int(p)}")

solve_price()
>>>4