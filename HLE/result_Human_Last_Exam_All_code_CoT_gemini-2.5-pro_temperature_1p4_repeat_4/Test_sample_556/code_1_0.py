import pandas as pd

def solve_price():
    """
    Solves for the price p based on the given principal-agent model and data.
    """
    # Step 1: Define the given data
    data = {
        'e': [22, 22, 10, 10],
        'y': [132, 44, 60, 20],
        'prob': [0.4375, 0.0625, 0.0625, 0.4375]
    }
    df = pd.DataFrame(data)

    # Step 2: Identify the states of the world, s_H and s_L
    # Calculate s = y / e
    df['s'] = df['y'] / df['e']
    s_values = sorted(df['s'].unique())
    s_L, s_H = s_values[0], s_values[1]

    e_values = sorted(df['e'].unique(), reverse=True)
    e_H, e_L = e_values[0], e_values[1]
    
    # Step 3: Map (e,y) probabilities to (theta,s) probabilities
    # We assume e_H (22) is chosen on signal theta_H and e_L (10) on signal theta_L
    P_thetaH_sH = df[(df['e'] == e_H) & (df['s'] == s_H)]['prob'].iloc[0]
    P_thetaH_sL = df[(df['e'] == e_H) & (df['s'] == s_L)]['prob'].iloc[0]
    P_thetaL_sH = df[(df['e'] == e_L) & (df['s'] == s_H)]['prob'].iloc[0]
    P_thetaL_sL = df[(df['e'] == e_L) & (df['s'] == s_L)]['prob'].iloc[0]
    
    # Step 4: Calculate marginal probabilities of signals
    P_thetaH = P_thetaH_sH + P_thetaH_sL
    P_thetaL = P_thetaL_sH + P_thetaL_sL

    # Step 5 & 6: Calculate conditional expectations E[s|theta]
    # E[s|theta_H] = s_H * P(s_H|theta_H) + s_L * P(s_L|theta_H)
    # P(s_H|theta_H) = P(theta_H, s_H) / P(theta_H)
    E_s_given_thetaH = (s_H * P_thetaH_sH + s_L * P_thetaH_sL) / P_thetaH
    
    # E[s|theta_L] = s_H * P(s_H|theta_L) + s_L * P(s_L|theta_L)
    # P(s_H|theta_L) = P(theta_L, s_H) / P(theta_L)
    E_s_given_thetaL = (s_H * P_thetaL_sH + s_L * P_thetaL_sL) / P_thetaL
    
    # Step 7: Use the employee's effort rule e = beta * E[s|theta] to find beta
    # We can use either the high effort or low effort rule.
    beta = e_H / E_s_given_thetaH
    
    # We can verify with the other rule: beta_check = e_L / E_s_given_thetaL
    # print(f"Beta check: {e_L} / {E_s_given_thetaL} = {beta_check}")
    
    # Step 8: The firm's profit maximization implies p = beta
    p = beta
    
    # Final Output
    print("The problem is solved by deriving the relationship p = β from the firm's profit maximization problem.")
    print("β is the slope of the wage contract, w = α + βy.")
    print("The employee's optimal effort is e = β * E[s|θ]. We can find β from the data.")
    print("\nCalculations:")
    print(f"High state s_H = {s_H}, Low state s_L = {s_L}")
    print(f"High effort e_H = {e_H}, Low effort e_L = {e_L}")
    print(f"Expected state given high signal, E[s|θ_H] = {E_s_given_thetaH}")
    print(f"The equation for β is: β = e_H / E[s|θ_H]")
    print(f"β = {e_H} / {E_s_given_thetaH} = {beta}")
    print(f"\nThe firm's optimal pricing leads to the condition p = β.")
    print(f"Therefore, the final equation for p is:")
    print(f"p = {p}")

solve_price()