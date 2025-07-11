import pandas as pd

def solve_for_p():
    """
    This function calculates the price p based on the provided probability distribution
    of effort and output.
    """
    # Step 1: Define the given probability distribution P(e,y) and determine the states 's'.
    # P(e,y) is the probability of observing effort 'e' and output 'y'.
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # The state 's' is calculated as y/e.
    # We build a DataFrame to facilitate calculations.
    data = []
    for (e, y), prob in prob_dist.items():
        s = y / e
        data.append({'e': e, 'y': y, 's': s, 'prob': prob})
    df = pd.DataFrame(data)

    # The two possible states and efforts
    s_values = sorted(df['s'].unique())
    s_L, s_H = s_values[0], s_values[1]
    e_values = sorted(df['e'].unique())
    e_L, e_H = e_values[0], e_values[1]

    print(f"From the data, we deduce the following model parameters:")
    print(f"High state s_H = {s_H}, Low state s_L = {s_L}")
    print(f"High effort e_H = {e_H}, Low effort e_L = {e_L}")
    print("-" * 30)

    # Step 2 & 3: Calculate conditional expectations E[s|e]
    # This corresponds to E[s|signal] that leads to effort e.

    # Probability of high effort e_H
    P_eH = df[df['e'] == e_H]['prob'].sum()
    # Probability of low effort e_L
    P_eL = df[df['e'] == e_L]['prob'].sum()

    # Calculate E[s | e_H]
    P_sH_given_eH = df[(df['e'] == e_H) & (df['s'] == s_H)]['prob'].sum() / P_eH
    P_sL_given_eH = df[(df['e'] == e_H) & (df['s'] == s_L)]['prob'].sum() / P_eH
    E_s_given_eH = s_H * P_sH_given_eH + s_L * P_sL_given_eH

    # Calculate E[s | e_L]
    P_sH_given_eL = df[(df['e'] == e_L) & (df['s'] == s_H)]['prob'].sum() / P_eL
    P_sL_given_eL = df[(df['e'] == e_L) & (df['s'] == s_L)]['prob'].sum() / P_eL
    E_s_given_eL = s_H * P_sH_given_eL + s_L * P_sL_given_eL
    
    print("The employee's optimal effort 'e' is given by e = \u03B2 \u00B7 E[s|signal].")
    print(f"For high effort: E[s|e_H] = {s_H}*{P_sH_given_eH:.3f} + {s_L}*{P_sL_given_eH:.3f} = {E_s_given_eH}")
    print(f"For low effort: E[s|e_L] = {s_H}*{P_sH_given_eL:.3f} + {s_L}*{P_sL_given_eL:.3f} = {E_s_given_eL}")
    print("-" * 30)

    # Step 4: Solve for beta (β)
    # The equation e = β * E[s|e] gives us beta.
    # We can use either the high-effort or low-effort case. Both should be consistent.
    beta = e_H / E_s_given_eH
    
    print("We can now solve for the contract parameter \u03B2 (beta).")
    print(f"Using the high effort case: \u03B2 = e_H / E[s|e_H]")
    print(f"\u03B2 = {e_H} / {E_s_given_eH} = {beta}")
    # Verifying with the low effort case:
    # beta_check = e_L / E_s_given_eL
    # print(f"Check with low effort case: {e_L} / {E_s_given_eL} = {beta_check}")
    print("-" * 30)
    
    # Step 5: State the economic result p = β
    # The firm maximizes profit E[p*y - w]. For a risk-neutral agent, the optimal
    # piece rate β equals the price p.
    p = beta
    
    print("To maximize profit, the firm sets the piece-rate \u03B2 equal to the output price p.")
    print("Therefore, p = \u03B2.")
    print(f"The final equation for p is:")
    print(f"p = \u03B2 = {e_H} / {E_s_given_eH} = {p}")


solve_for_p()