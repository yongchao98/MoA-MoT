import collections

def solve_for_p():
    """
    Calculates the value of p based on the provided probability distribution.
    """
    # 1. Deconstruct the Data
    # The given probability distribution of pairs (e, y)
    P_ey = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375,
    }

    print("Step 1: Determine states (s) and effort levels (e) from the data.\n")
    
    # Extract unique effort and state values
    efforts = sorted(list(set(e for e, y in P_ey.keys())), reverse=True)
    e_H, e_L = efforts[0], efforts[1]
    
    states = sorted(list(set(y / e for e, y in P_ey.keys())), reverse=True)
    s_H, s_L = states[0], states[1]

    print(f"The high effort level e_H is {e_H} and the low effort level e_L is {e_L}.")
    print(f"The high state of the world s_H is {s_H} and the low state s_L is {s_L}.\n")

    # 2. Calculate Conditional Expectations
    print("Step 2: Calculate the employee's conditional expectations of the state, E[s|θ].\n")

    # Map (e,y) probabilities to (s,θ) probabilities
    # Signal s_H (θ_H) leads to effort e_H. Signal s_L (θ_L) leads to effort e_L.
    # P(s=s_H, θ=s_H) = P(e=e_H, y=s_H*e_H)
    P_sH_thetaH = P_ey[(e_H, s_H * e_H)]
    # P(s=s_L, θ=s_H) = P(e=e_H, y=s_L*e_H)
    P_sL_thetaH = P_ey[(e_H, s_L * e_H)]
    # P(s=s_H, θ=s_L) = P(e=e_L, y=s_H*e_L)
    P_sH_thetaL = P_ey[(e_L, s_H * e_L)]
    # P(s=s_L, θ=s_L) = P(e=e_L, y=s_L*e_L)
    P_sL_thetaL = P_ey[(e_L, s_L * e_L)]

    # Calculate marginal probabilities of signals
    P_thetaH = P_sH_thetaH + P_sL_thetaH
    P_thetaL = P_sH_thetaL + P_sL_thetaL

    # Calculate conditional probabilities P(s|θ) using Bayes' rule
    P_sH_given_thetaH = P_sH_thetaH / P_thetaH
    P_sL_given_thetaH = P_sL_thetaH / P_thetaH
    P_sH_given_thetaL = P_sH_thetaL / P_thetaL
    P_sL_given_thetaL = P_sL_thetaL / P_thetaL

    # Calculate conditional expectations E[s|θ]
    E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
    E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL

    print("The probability of receiving signal θ_H is P(θ=s_H) = "
          f"{P_sH_thetaH} + {P_sL_thetaH} = {P_thetaH}")
    print("The conditional expectation of s given signal θ_H is:")
    print(f"E[s|θ=s_H] = s_H * P(s=s_H|θ=s_H) + s_L * P(s=s_L|θ=s_H)")
    print(f"E[s|θ=s_H] = {s_H} * ({P_sH_thetaH}/{P_thetaH}) + {s_L} * ({P_sL_thetaH}/{P_thetaH})")
    print(f"E[s|θ=s_H] = {s_H} * {P_sH_given_thetaH} + {s_L} * {P_sL_given_thetaH} = {E_s_given_thetaH}\n")

    print("The probability of receiving signal θ_L is P(θ=s_L) = "
          f"{P_sH_thetaL} + {P_sL_thetaL} = {P_thetaL}")
    print("The conditional expectation of s given signal θ_L is:")
    print(f"E[s|θ=s_L] = s_H * P(s=s_H|θ=s_L) + s_L * P(s=s_L|θ=s_L)")
    print(f"E[s|θ=s_L] = {s_H} * ({P_sH_thetaL}/{P_thetaL}) + {s_L} * ({P_sL_thetaL}/{P_thetaL})")
    print(f"E[s|θ=s_L] = {s_H} * {P_sH_given_thetaL} + {s_L} * {P_sL_given_thetaL} = {E_s_given_thetaL}\n")

    # 3. Determine β
    print("Step 3: Determine the contract parameter β.\n")
    print("The employee's optimal effort is e(θ) = β * E[s|θ]. We can find β using the high effort case.")
    # beta = e_H / E_s_given_thetaH
    beta = e_H / E_s_given_thetaH
    print(f"β = e_H / E[s|θ=s_H] = {e_H} / {E_s_given_thetaH}")
    print(f"β = {beta}\n")

    # (Optional) Check for consistency
    beta_check = e_L / E_s_given_thetaL
    print(f"Checking with the low effort case: β = e_L / E[s|θ=s_L] = {e_L} / {E_s_given_thetaL} = {beta_check}. The values match.\n")
    
    # 4. Determine p
    print("Step 4: Relate β to the price p.\n")
    print("In a principal-agent model with a risk-neutral agent, the profit-maximizing firm sets the incentive parameter β equal to the output price p.")
    print("Therefore, p = β.\n")

    # 5. Final Answer
    p = beta
    print("Final Answer:")
    print(f"The value of p is equal to the calculated value of β, which is {p}.")

solve_for_p()
<<<4>>>