def solve_economy_shock():
    """
    Calculates the percentage change in nominal wage, price of good X,
    and consumption of good X based on the given economic model and shock.
    """
    
    # --- Given Parameters ---
    # Cost shares
    theta_LX = 2/3  # Labor's share in X
    theta_KX = 1/3  # Capital's share in X
    theta_LY = 1/3  # Labor's share in Y
    theta_KY = 2/3  # Capital's share in Y

    # Demand elasticities
    eta_X = 1.0     # Income elasticity for X
    eta_Y = 1.0     # Income elasticity for Y
    epsilon_XX = -2.0 # Own-price elasticity for X
    epsilon_YY = -1.0 # Own-price elasticity for Y
    
    # Shocks and known outcomes
    tax_rate = 0.02       # 2% tax, so r_X changes by 0.02
    C_Y_hat = 0.03        # Consumption of Y increases by 3%
    
    # --- Step-by-step Calculation ---
    
    # 1. Calculate percentage change in nominal wage (w_hat)
    # The price of good Y (P_Y) and return to capital for sector Y (r_Y) are fixed
    # by the world market.
    # P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat
    # 0 = (1/3) * w_hat + (2/3) * 0
    P_Y_hat = 0.0
    r_Y_hat = 0.0
    w_hat = (P_Y_hat - theta_KY * r_Y_hat) / theta_LY
    
    # 2. Calculate percentage change in price of good X (P_X_hat)
    # The tax increases the cost of capital for producers in sector X.
    # P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
    r_X_hat = tax_rate
    P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat

    # 3. Calculate percentage change in income (I_hat)
    # First, find cross-price elasticity for Y using homogeneity of demand.
    # epsilon_YX + epsilon_YY + eta_Y = 0
    epsilon_YX = -(epsilon_YY + eta_Y)
    
    # Use the demand equation for Y to find income change.
    # C_Y_hat = epsilon_YX * P_X_hat + epsilon_YY * P_Y_hat + eta_Y * I_hat
    I_hat = (C_Y_hat - epsilon_YX * P_X_hat - epsilon_YY * P_Y_hat) / eta_Y

    # 4. Calculate percentage change in consumption of good X (C_X_hat)
    # First, find cross-price elasticity for X using homogeneity of demand.
    # epsilon_XX + epsilon_XY + eta_X = 0
    epsilon_XY = -(epsilon_XX + eta_X)
    
    # Use the demand equation for X.
    # C_X_hat = epsilon_XX * P_X_hat + epsilon_XY * P_Y_hat + eta_X * I_hat
    C_X_hat = epsilon_XX * P_X_hat + epsilon_XY * P_Y_hat + eta_X * I_hat

    # Convert results to percentages for output
    w_hat_percent = w_hat * 100
    P_X_hat_percent = P_X_hat * 100
    C_X_hat_percent = C_X_hat * 100

    # Print the final result as 3 comma-separated values
    print(f"{w_hat_percent},{P_X_hat_percent},{C_X_hat_percent}")

solve_economy_shock()