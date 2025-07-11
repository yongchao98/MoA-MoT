import collections

def solve_for_p():
    """
    This function solves for the price 'p' based on the provided economic model and data.
    """
    
    # 1. Define the given probabilities and observed values
    # P(e, y)
    probs = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375,
    }

    # Observed efforts. The higher effort corresponds to the high signal (theta_H)
    e_H = 22
    e_L = 10

    print("Step 1: Determine the states of the world (s_H, s_L) from the data.")
    # For e=22, y can be 132 or 44. Since y = s*e, s = y/e.
    s_values = {y / e for e, y in probs.keys() if e == e_H}
    s_H = max(s_values)
    s_L = min(s_values)
    print(f"The high state is s_H = 132 / 22 = {s_H}")
    print(f"The low state is s_L = 44 / 22 = {s_L}")
    print("-" * 30)

    # 2. Calculate probabilities related to signals and states.
    # The event (e,y) corresponds to a unique (signal, state) pair.
    # e=e_H (22) implies signal theta_H. e=e_L (10) implies signal theta_L.
    # s=s_H (6) implies state s_H. s=s_L (2) implies state s_L.
    
    # Joint probabilities of signal and state
    P_thetaH_sH = probs[(e_H, s_H * e_H)]  # P(signal=H, state=H)
    P_thetaH_sL = probs[(e_H, s_L * e_H)]  # P(signal=H, state=L)
    P_thetaL_sH = probs[(e_L, s_H * e_L)]  # P(signal=L, state=H)
    P_thetaL_sL = probs[(e_L, s_L * e_L)]  # P(signal=L, state=L)

    # Marginal probability of each signal
    P_thetaH = P_thetaH_sH + P_thetaH_sL
    P_thetaL = P_thetaL_sH + P_thetaL_sL

    print("Step 2: Calculate the employee's conditional expectations of the state.")
    print(f"Probability of receiving high signal P(theta=H) = {P_thetaH_sH} + {P_thetaH_sL} = {P_thetaH}")
    print(f"Probability of receiving low signal P(theta=L) = {P_thetaL_sH} + {P_thetaL_sL} = {P_thetaL}")

    # Conditional probabilities of state given signal
    P_sH_given_thetaH = P_thetaH_sH / P_thetaH
    P_sL_given_thetaH = P_thetaH_sL / P_thetaH
    
    P_sH_given_thetaL = P_thetaL_sH / P_thetaL
    P_sL_given_thetaL = P_thetaL_sL / P_thetaL

    # Employee's conditional expectations of the state
    E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
    E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL
    
    print(f"Employee's expectation of s given a high signal, E[s|theta_H] = {s_H}*{P_sH_given_thetaH:.3f} + {s_L}*{P_sL_given_thetaH:.3f} = {E_s_given_thetaH}")
    print(f"Employee's expectation of s given a low signal, E[s|theta_L] = {s_H}*{P_sH_given_thetaL:.3f} + {s_L}*{P_sL_given_thetaL:.3f} = {E_s_given_thetaL}")
    print("-" * 30)

    # 3. Use the employee's optimal effort rule to find beta.
    # The employee's optimal effort is e_theta = beta * E[s|theta].
    # We can use the high signal case to find beta.
    beta = e_H / E_s_given_thetaH
    
    # We can check this with the low signal case as well.
    beta_check = e_L / E_s_given_thetaL

    print("Step 3: Determine the contract parameter beta.")
    print("The employee chooses effort to maximize E[alpha + beta*s*e - e^2/2 | theta], which gives e = beta * E[s|theta].")
    print(f"Using the high signal case: beta = e_H / E[s|theta_H]")
    print(f"The equation is: beta = {e_H} / {E_s_given_thetaH}")
    print(f"So, beta = {beta}")
    print(f"(Check with low signal case: beta = {e_L} / {E_s_given_thetaL} = {beta_check})")
    print("-" * 30)

    # 4. Use the firm's optimality condition to find p.
    # The firm chooses beta to maximize total surplus E[p*y - e^2/2].
    # This optimization leads to the condition p = beta.
    p = beta

    print("Step 4: Determine the price p.")
    print("The firm chooses beta to maximize expected profit. In this model, this implies maximizing total surplus.")
    print("The optimal choice for the firm leads to the condition: p = beta.")
    print(f"Therefore, the value of p is {p}.")

if __name__ == '__main__':
    solve_for_p()
    # Final Answer
    # The code calculates that beta is 4. Since p = beta, p = 4.
    p_final = 4
    print(f"\n<<<_p_>>>{p_final}<<<p>>>")