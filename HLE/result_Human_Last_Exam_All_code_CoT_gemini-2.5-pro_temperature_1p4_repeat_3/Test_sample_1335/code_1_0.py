def solve_economy_tax_effect():
    """
    Calculates the percentage change in nominal wage, price of good X, and consumption of good X
    based on the given economic model and a tax on capital in sector X.
    """

    # Given parameters
    hat_r_X = 0.02  # 2% tax on capital return in sector X, a 2% increase in firm cost
    hat_C_Y = 0.03  # Consumption of Y increases by 3%

    # Factor shares of cost
    theta_LX = 2/3
    theta_KX = 1/3
    theta_LY = 1/3
    theta_KY = 2/3

    # Demand elasticities
    eta_X = 1.0
    eta_Y = 1.0
    eps_XX = -2.0
    eps_YY = -1.0
    
    # --- Step 1: Solve for percentage change in wage (hat_w) ---
    # From the zero-profit condition in the traded good sector Y:
    # hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y
    # As a small open economy, hat_P_Y = 0 and hat_r_Y = 0.
    hat_w = 0.0

    # --- Step 2: Solve for percentage change in price of X (hat_P_X) ---
    # From the zero-profit condition in sector X:
    # hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
    hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
    
    # --- Step 3: Solve for percentage change in nominal income (hat_I) ---
    # From the demand function for good Y:
    # hat_C_Y = epsilon_YX * hat_P_X + epsilon_YY * hat_P_Y + eta_Y * hat_I
    # The given elasticities for Y (eps_YY=-1, eta_Y=1) imply a cross-price 
    # elasticity epsilon_YX = 0 to satisfy the homogeneity condition.
    # Also, hat_P_Y = 0.
    # So, hat_C_Y = eta_Y * hat_I
    hat_I = hat_C_Y / eta_Y

    # --- Step 4: Solve for percentage change in consumption of X (hat_C_X) ---
    # From the demand function for good X:
    # hat_C_X = eps_XX * hat_P_X + epsilon_XY * hat_P_Y + eta_X * hat_I
    # With hat_P_Y = 0, the equation simplifies.
    hat_C_X = eps_XX * hat_P_X + eta_X * hat_I

    # --- Step 5: Convert results to percentages and print ---
    # The question asks for the percentage change on nominal wage, price of good X, and consumption of good X.
    # Multiplying the proportional changes by 100 gives the percentage values.
    
    wage_change_percent = hat_w * 100
    price_X_change_percent = hat_P_X * 100
    consumption_X_change_percent = hat_C_X * 100

    print(f"The tax on capital in sector X is {hat_r_X*100}%.")
    print(f"As a result, the consumption of Y increases by {hat_C_Y*100}%.")
    print("Based on these changes, we can calculate the following effects:")
    print("---")
    print(f"1. The percentage change in the nominal wage is:")
    print(f"{hat_w} * 100 = {wage_change_percent}%")
    print("---")
    print(f"2. The percentage change in the price of good X is:")
    print(f"({theta_LX:.2f} * {hat_w} + {theta_KX:.2f} * {hat_r_X}) * 100 = {price_X_change_percent:.4f}%")
    print("---")
    print(f"3. The percentage change in the consumption of good X is:")
    print(f"({eps_XX} * {hat_P_X:.4f} + {eta_X} * {hat_I}) * 100 = {consumption_X_change_percent:.4f}%")
    print("---")
    print("Final comma-separated values (Wage %, Px %, Cx %):")
    print(f"{wage_change_percent},{price_X_change_percent:.4f},{consumption_X_change_percent:.4f}")

solve_economy_tax_effect()
<<<0.0,0.6667,1.6667>>>