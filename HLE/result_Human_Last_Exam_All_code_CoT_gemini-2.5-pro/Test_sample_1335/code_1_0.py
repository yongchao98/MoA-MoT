def solve_economy_tax_effect():
    """
    Calculates the percentage change on nominal wage, price of good X, 
    and consumption of good X due to a tax on capital in sector X.
    """
    # --- Given Parameters ---
    # Cost shares for labor (L) and capital (K) in sectors X and Y
    theta_LX = 2/3
    theta_KX = 1/3
    theta_LY = 1/3
    theta_KY = 2/3

    # Own-price and income elasticities of demand
    eta_XX = -2
    eta_IX = 1
    eta_IY = 1

    # The shock: a 2% tax on capital return in sector X
    r_X_hat = 0.02

    # The observed outcome: consumption of Y increases by 3%
    C_Y_hat = 0.03
    
    # In a small open economy, the price of traded goods and the return to mobile
    # factors are fixed by the world market.
    P_Y_hat = 0
    r_Y_hat = 0

    # --- Step 1: Calculate the percentage change in nominal wage (w_hat) ---
    # From the zero-profit condition in sector Y: P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat
    # 0 = (1/3) * w_hat + (2/3) * 0
    w_hat = 0

    # --- Step 2: Calculate the percentage change in the price of good X (P_X_hat) ---
    # From the zero-profit condition in sector X: P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
    P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat

    # --- Step 3: Calculate the percentage change in nominal income (I_hat) ---
    # From the demand function for Y: C_Y_hat = eta_YX * P_X_hat + eta_YY * P_Y_hat + eta_IY * I_hat
    # We assume the cross-price elasticity (eta_YX) is 0 as it's not provided.
    # 0.03 = 0 * P_X_hat + (-1) * 0 + 1 * I_hat
    I_hat = C_Y_hat / eta_IY

    # --- Step 4: Calculate the percentage change in consumption of good X (C_X_hat) ---
    # From the demand function for X: C_X_hat = eta_XX * P_X_hat + eta_XY * P_Y_hat + eta_IX * I_hat
    # We assume cross-price elasticity (eta_XY) is 0.
    C_X_hat = eta_XX * P_X_hat + eta_IX * I_hat

    # Convert results to percentages
    w_hat_pct = w_hat * 100
    P_X_hat_pct = P_X_hat * 100
    C_X_hat_pct = C_X_hat * 100

    # Print the final comma-separated values
    # The final equation for the wage is 0 = (1/3) * w_hat, so w_hat = 0%
    # The final equation for the price of X is P_X_hat = (1/3) * 2% = 2/3%
    # The final equation for consumption of X is C_X_hat = -2 * (2/3%) + 1 * 3% = -4/3% + 9/3% = 5/3%
    print(f"Percentage change in nominal wage: {w_hat_pct}")
    print(f"Percentage change in price of good X: {P_X_hat_pct}")
    print(f"Percentage change in consumption of good X: {C_X_hat_pct}")
    print("\nFinal Answer:")
    print(f"{w_hat_pct},{P_X_hat_pct},{C_X_hat_pct}")

solve_economy_tax_effect()
<<<0,0.6666666666666666,1.6666666666666667>>>