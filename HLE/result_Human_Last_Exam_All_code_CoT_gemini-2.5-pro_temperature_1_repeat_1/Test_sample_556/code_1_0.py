import collections

def solve_for_p():
    """
    Calculates the value of p based on the given probability distribution.
    """
    # Step 1: Infer model parameters from the data
    prob_data = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # Efforts are e_H = 22 and e_L = 10
    e_H = 22
    e_L = 10

    # From y = s * e, we can find the states s.
    # For e = 22, y can be 132 or 44. s = 132/22=6 or s = 44/22=2.
    # So, s_H = 6 and s_L = 2.
    s_H = 132 / 22
    s_L = 44 / 22
    print(f"Inferred model parameters:")
    print(f"High effort e_H = {e_H}")
    print(f"Low effort e_L = {e_L}")
    print(f"High state s_H = {s_H}")
    print(f"Low state s_L = {s_L}\n")

    # Step 2: Analyze Employee's Behavior to find beta
    # The four outcomes correspond to (signal, state) pairs:
    # (e=22, y=132) -> (signal=s_H, state=s_H)
    # (e=22, y=44)  -> (signal=s_H, state=s_L)
    # (e=10, y=60)  -> (signal=s_L, state=s_H)
    # (e=10, y=20)  -> (signal=s_L, state=s_L)

    # Calculate marginal probability of signals
    prob_signal_H = prob_data[(22, 132)] + prob_data[(22, 44)]
    prob_signal_L = prob_data[(10, 60)] + prob_data[(10, 20)]

    # Calculate conditional probabilities of states given signals
    prob_sH_given_signalH = prob_data[(22, 132)] / prob_signal_H
    prob_sL_given_signalH = prob_data[(22, 44)] / prob_signal_H
    
    prob_sH_given_signalL = prob_data[(10, 60)] / prob_signal_L
    prob_sL_given_signalL = prob_data[(10, 20)] / prob_signal_L

    # Calculate conditional expectations E[s | signal]
    E_s_given_signalH = s_H * prob_sH_given_signalH + s_L * prob_sL_given_signalH
    E_s_given_signalL = s_H * prob_sH_given_signalL + s_L * prob_sL_given_signalL

    print("Calculating employee's conditional expectations:")
    print(f"E[s | signal=s_H] = {s_H} * {prob_sH_given_signalH:.4f} + {s_L} * {prob_sL_given_signalH:.4f} = {E_s_given_signalH}")
    print(f"E[s | signal=s_L] = {s_H} * {prob_sH_given_signalL:.4f} + {s_L} * {prob_sL_given_signalL:.4f} = {E_s_given_signalL}\n")

    # Use the employee's FOC: e = beta * E[s | signal] to find beta
    beta_from_H = e_H / E_s_given_signalH
    beta_from_L = e_L / E_s_given_signalL
    
    # The value of beta should be consistent
    beta = beta_from_H

    print("Calculating beta from employee's First-Order Condition (e = beta * E[s|signal]):")
    print(f"From high effort: beta = {e_H} / {E_s_given_signalH} = {beta}")
    print(f"From low effort:  beta = {e_L} / {E_s_given_signalL} = {beta_from_L}\n")
    
    # Step 3 & 4: Use the firm's optimization result p = beta to find p
    # The firm's profit maximization FOC for beta results in p = beta.
    p = beta
    
    print("The firm's profit maximization implies a relationship between price 'p' and contract parameter 'beta'.")
    print("For a linear contract and quadratic effort costs, this relationship is p = beta.")
    print("\nFinal calculation for p:")
    print(f"p = \u03B2 = {p}")

solve_for_p()