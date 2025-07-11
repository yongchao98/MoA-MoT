import pandas as pd

def solve_price():
    """
    Calculates the price 'p' based on the provided probability distribution
    of effort and output.
    """
    # 1. Define the given probability distribution
    data = {
        'e': [22, 22, 10, 10],
        'y': [132, 44, 60, 20],
        'prob': [0.4375, 0.0625, 0.0625, 0.4375]
    }
    df = pd.DataFrame(data)

    # 2. Determine the states of the world 's' using y = s * e
    df['s'] = df['y'] / df['e']
    s_values = sorted(df['s'].unique())
    s_H, s_L = s_values[1], s_values[0]

    print("Step 1: Determine the states of the world 's' and effort levels 'e'.")
    print(f"From the (e, y) pairs, we find two possible states of the world:")
    print(f"s_H = {s_H}")
    print(f"s_L = {s_L}")

    e_values = sorted(df['e'].unique(), reverse=True)
    e_H, e_L = e_values[0], e_values[1]
    print(f"\nThe two observed effort levels are:")
    print(f"e_H = {e_H} (effort after a high signal)")
    print(f"e_L = {e_L} (effort after a low signal)\n")
    print("---------------------\n")

    # 3. Calculate conditional expectations and beta
    print("Step 2: Use the employee's optimization rule e = β * E[s|θ] to find β.")
    
    # Probabilities for calculations
    p_eH = df[df['e'] == e_H]['prob'].sum()
    p_sH_given_eH = df[(df['e'] == e_H) & (df['s'] == s_H)]['prob'].iloc[0] / p_eH
    p_sL_given_eH = df[(df['e'] == e_H) & (df['s'] == s_L)]['prob'].iloc[0] / p_eH

    # E[s | θ_H]
    E_s_given_eH = s_H * p_sH_given_eH + s_L * p_sL_given_eH
    print(f"The conditional expectation of s given a high signal (e_H={e_H}) is:")
    print(f"E[s|θ_H] = {s_H}*P(s=s_H|θ_H) + {s_L}*P(s=s_L|θ_H) = {s_H}*{p_sH_given_eH:.4f} + {s_L}*{p_sL_given_eH:.4f} = {E_s_given_eH}")
    
    # Calculate beta from the high effort case
    beta_from_eH = e_H / E_s_given_eH
    print(f"From e_H = β * E[s|θ_H], we get:")
    print(f"{e_H} = β * {E_s_given_eH}")
    print(f"β = {e_H} / {E_s_given_eH} = {beta_from_eH:.1f}\n")

    # Repeat for low effort to confirm
    p_eL = df[df['e'] == e_L]['prob'].sum()
    p_sH_given_eL = df[(df['e'] == e_L) & (df['s'] == s_H)]['prob'].iloc[0] / p_eL
    p_sL_given_eL = df[(df['e'] == e_L) & (df['s'] == s_L)]['prob'].iloc[0] / p_eL
    
    # E[s | θ_L]
    E_s_given_eL = s_H * p_sH_given_eL + s_L * p_sL_given_eL
    print(f"To confirm, the conditional expectation of s given a low signal (e_L={e_L}) is:")
    print(f"E[s|θ_L] = {s_H}*P(s=s_H|θ_L) + {s_L}*P(s=s_L|θ_L) = {s_H}*{p_sH_given_eL:.4f} + {s_L}*{p_sL_given_eL:.4f} = {E_s_given_eL}")
    
    beta_from_eL = e_L / E_s_given_eL
    print(f"From e_L = β * E[s|θ_L], we get:")
    print(f"{e_L} = β * {E_s_given_eL}")
    print(f"β = {e_L} / {E_s_given_eL} = {beta_from_eL:.1f}")
    
    # Final beta value
    beta = beta_from_eH
    print(f"\nThe value of β is consistent across both effort levels. β = {beta:.1f}\n")
    print("---------------------\n")
    
    # 4. Use the firm's optimization to find p
    print("Step 3: Use the firm's optimization rule to find p.")
    print("The firm chooses β to maximize its expected profit. The first-order condition")
    print("of this problem simplifies to p = β.")
    
    p = beta
    print("\nFinal Equation and Answer:")
    print(f"The price p is equal to the contract parameter β.")
    print(f"p = β")
    print(f"Therefore, p = {p:.1f}")
    
    return p

# Execute the function to find the answer
p_value = solve_price()