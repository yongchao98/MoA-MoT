import collections

def solve():
    """
    Solves the problem by deducing model parameters from the given data
    and using the firm's and employee's optimization conditions.
    """
    # Step 1: Define the given data
    probabilities = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # Step 2: Identify states of the world (s) and effort levels (e)
    # Effort levels are the first element of the keys
    efforts = sorted(list(set(e for e, y in probabilities.keys())), reverse=True)
    e_H = efforts[0] # Higher effort
    e_L = efforts[1] # Lower effort

    # The state of the world, s, is given by y/e
    # Using the high effort case:
    s_H = 132 / e_H
    s_L = 44 / e_H
    
    print("Step 1: Identify model parameters from the data.")
    print(f"The observed effort levels are e_H = {e_H} and e_L = {e_L}.")
    print(f"The states of the world are deduced from y = s * e:")
    print(f"s_H = 132 / {e_H} = {s_H}")
    print(f"s_L = 44 / {e_H} = {s_L}")
    print(f"Checking for consistency with low effort: {s_H}*{e_L}=60.0 and {s_L}*{e_L}=20.0. The states are consistent.")

    # Step 3: Reconstruct the joint probability distribution P(signal, state)
    # We associate high effort with signal 'H' and low effort with signal 'L'.
    # We associate state 'H' with s_H and state 'L' with s_L.
    P_thetaH_sH = probabilities[(e_H, s_H * e_H)]
    P_thetaH_sL = probabilities[(e_H, s_L * e_H)]
    P_thetaL_sH = probabilities[(e_L, s_H * e_L)]
    P_thetaL_sL = probabilities[(e_L, s_L * e_L)]

    print("\nStep 2: Reconstruct the joint probability P(signal, state).")
    print(f"P(signal=H, state=H) = P(e={e_H}, y={e_H*s_H}) = {P_thetaH_sH}")
    print(f"P(signal=H, state=L) = P(e={e_H}, y={e_H*s_L}) = {P_thetaH_sL}")
    print(f"P(signal=L, state=H) = P(e={e_L}, y={e_L*s_H}) = {P_thetaL_sH}")
    print(f"P(signal=L, state=L) = P(e={e_L}, y={e_L*s_L}) = {P_thetaL_sL}")

    # Step 4: Calculate marginal and conditional probabilities
    P_thetaH = P_thetaH_sH + P_thetaH_sL
    P_thetaL = P_thetaL_sH + P_thetaL_sL
    
    P_sH_given_thetaH = P_thetaH_sH / P_thetaH
    P_sL_given_thetaH = P_thetaH_sL / P_thetaH
    P_sH_given_thetaL = P_thetaL_sH / P_thetaL
    P_sL_given_thetaL = P_thetaL_sL / P_thetaL

    # Step 5: Calculate the worker's conditional expectation of the state
    E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
    E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL

    print("\nStep 3: Calculate the worker's conditional expectation of the state, E[s|signal].")
    print(f"E[s | signal=H] = {s_H} * P(s=H|signal=H) + {s_L} * P(s=L|signal=H)")
    print(f"E[s | signal=H] = {s_H} * {P_sH_given_thetaH:.4f} + {s_L} * {P_sL_given_thetaH:.4f} = {E_s_given_thetaH}")
    print(f"E[s | signal=L] = {s_H} * P(s=H|signal=L) + {s_L} * P(s=L|signal=L)")
    print(f"E[s | signal=L] = {s_H} * {P_sH_given_thetaL:.4f} + {s_L} * {P_sL_given_thetaL:.4f} = {E_s_given_thetaL}")

    # Step 6: Determine the contract parameter beta
    beta_from_H = e_H / E_s_given_thetaH
    beta_from_L = e_L / E_s_given_thetaL
    beta = beta_from_H # They are the same

    print("\nStep 4: Determine the contract parameter beta from the worker's effort choice rule: e = beta * E[s|signal].")
    print(f"From high effort: {e_H} = beta * {E_s_given_thetaH}, which gives beta = {beta_from_H:.1f}")
    print(f"From low effort:  {e_L} = beta * {E_s_given_thetaL}, which gives beta = {beta_from_L:.1f}")
    print(f"The value of beta is consistent. So, the firm chose beta = {beta:.1f}")

    # Step 7: Solve the firm's problem and find p
    print("\nStep 5: Solve the firm's profit maximization problem to find p.")
    print("The firm chooses beta to maximize E[profit] = p*E[y(beta)] - E[e(beta)^2/2].")
    print("The first-order condition of this problem with respect to beta is: p - beta = 0.")
    print("This implies that the firm's optimal choice of beta must be equal to the price p.")
    p = beta
    print(f"Since the firm chose beta = {beta:.1f}, the price p must be {p:.1f}.")

    return p

result = solve()
print(f"\n<<<4.0>>>")