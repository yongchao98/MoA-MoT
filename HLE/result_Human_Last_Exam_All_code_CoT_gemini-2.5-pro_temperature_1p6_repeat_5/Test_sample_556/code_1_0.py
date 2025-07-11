import pandas as pd

def solve_price_p():
    """
    Solves for the price p based on the given economic model and data.
    """
    # Step 1: Parse the data and find state values
    # Given probability distribution P(e, y)
    P_ey = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # Identify unique effort levels
    efforts = sorted(list(set(e for e, y in P_ey.keys())), reverse=True)
    e_H, e_L = efforts[0], efforts[1]

    # For e_H, find the possible values of s = y/e
    y_for_eH = [y for e, y in P_ey.keys() if e == e_H]
    s_values = sorted([y / e_H for y in y_for_eH], reverse=True)
    s_H, s_L = s_values[0], s_values[1]

    print(f"Step 1: Analyzing the data.")
    print(f"The high effort level is e_H = {e_H} and the low effort level is e_L = {e_L}.")
    print(f"From the relationship y = s * e, we find the states of the world: s_H = {s_H} and s_L = {s_L}.")
    print("-" * 30)

    # Step 2: Calculate worker's conditional expectations
    # P(e=e_H) is the sum of probabilities for all outcomes with e=e_H
    P_e_H = P_ey.get((e_H, s_H * e_H), 0) + P_ey.get((e_H, s_L * e_H), 0)
    # P(e=e_L) is the sum of probabilities for all outcomes with e=e_L
    P_e_L = P_ey.get((e_L, s_H * e_L), 0) + P_ey.get((e_L, s_L * e_L), 0)

    # P(s=s_H | e=e_H) = P(e=e_H, s=s_H) / P(e=e_H)
    P_sH_given_eH = P_ey[(e_H, s_H * e_H)] / P_e_H
    P_sL_given_eH = P_ey[(e_H, s_L * e_H)] / P_e_H
    
    # E[s | signal leading to e_H]
    E_s_given_eH = s_H * P_sH_given_eH + s_L * P_sL_given_eH

    # P(s=s_H | e=e_L) = P(e=e_L, s=s_H) / P(e=e_L)
    P_sH_given_eL = P_ey[(e_L, s_H * e_L)] / P_e_L
    P_sL_given_eL = P_ey[(e_L, s_L * e_L)] / P_e_L

    # E[s | signal leading to e_L]
    E_s_given_eL = s_H * P_sH_given_eL + s_L * P_sL_given_eL
    
    print("Step 2: Calculating worker's beliefs (conditional expectations).")
    print(f"The worker's expected value of s after a high signal (leading to e={e_H}) is: E[s|θ_H] = {s_H}*P(s_H|θ_H) + {s_L}*P(s_L|θ_H) = {E_s_given_eH}")
    print(f"The worker's expected value of s after a low signal (leading to e={e_L}) is: E[s|θ_L] = {s_H}*P(s_H|θ_L) + {s_L}*P(s_L|θ_L) = {E_s_given_eL}")
    print("-" * 30)
    
    # Step 3: Determine the contract parameter beta
    # From the worker's effort choice e = beta * E[s|signal]
    beta_from_eH = e_H / E_s_given_eH
    beta_from_eL = e_L / E_s_given_eL
    
    # The value of beta should be consistent
    beta = beta_from_eH
    
    print("Step 3: Determining the contract parameter β.")
    print("From the worker's optimal effort rule: e = β * E[s|signal].")
    print(f"Using the high signal case: {e_H} = β * {E_s_given_eH}  =>  β = {beta_from_eH}")
    print(f"Using the low signal case: {e_L} = β * {E_s_given_eL}  =>  β = {beta_from_eL}")
    print(f"The data is consistent with β = {beta}.")
    print("-" * 30)

    # Step 4: Analyze the firm's problem and solve for p
    # The firm's profit maximization problem with respect to beta leads to the
    # first-order condition: p = beta.
    p = beta
    
    print("Step 4: Solving for the price p.")
    print("The firm chooses β to maximize its expected profit. The first-order condition of this problem simplifies to p = β.")
    print(f"The final equation is p = β.")
    print(f"Since β = {beta}, the price p must be {p}.")
    print("-" * 30)
    print(f"The final answer is p = {p}")

solve_price_p()
<<<4>>>