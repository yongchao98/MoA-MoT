def solve_economy_shock():
    """
    This script calculates the percentage change in nominal wage, the price of good X,
    and the consumption of good X in response to a tax on capital in sector X.
    """

    # --- Step 1: Define given parameters ---
    # The shock is a 2% tax on the return to capital in sector X.
    # This means the percentage change in the cost of capital for producers in sector X is 2%.
    r_X_hat = 0.02

    # The resulting change in the consumption of good Y is a 3% increase.
    d_Y_hat = 0.03

    # Cost shares
    theta_LX = 2/3  # Labor's share in X
    theta_KX = 1/3  # Capital's share in X
    theta_LY = 1/3  # Labor's share in Y
    theta_KY = 2/3  # Capital's share in Y

    # Demand elasticities
    eta_I_X = 1      # Income elasticity for X
    eta_I_Y = 1      # Income elasticity for Y
    epsilon_XX = -2  # Own-price elasticity for X
    epsilon_YY = -1  # Own-price elasticity for Y

    # --- Step 2: Use model assumptions to find wage change (w_hat) ---
    # In a small open economy, the price of the traded good Y is fixed (P_Y_hat = 0).
    P_Y_hat = 0
    # With international capital mobility, the world return to capital is fixed (r_hat = 0).
    # The tax is only in sector X, so the cost of capital in Y is unchanged (r_Y_hat = 0).
    r_Y_hat = 0

    # The zero-profit condition in sector Y is: P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat
    # 0 = (1/3) * w_hat + (2/3) * 0
    w_hat = 0
    print(f"The percentage change in nominal wage (w_hat) is derived from the zero-profit condition in sector Y:")
    print(f"P_Y_hat = θ_LY * w_hat + θ_KY * r_Y_hat")
    print(f"{P_Y_hat} = {theta_LY:.2f} * w_hat + {theta_KY:.2f} * {r_Y_hat}")
    print(f"Solving for w_hat gives: {w_hat * 100:.2f}%\n")


    # --- Step 3: Find price change of good X (P_X_hat) ---
    # The zero-profit condition in sector X is: P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
    P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
    print(f"The percentage change in the price of good X (P_X_hat) is derived from the zero-profit condition in sector X:")
    print(f"P_X_hat = θ_LX * w_hat + θ_KX * r_X_hat")
    print(f"P_X_hat = {theta_LX:.2f} * {w_hat * 100:.2f}% + {theta_KX:.2f} * {r_X_hat * 100:.2f}%")
    print(f"Solving for P_X_hat gives: {P_X_hat * 100:.2f}%\n")

    # --- Step 4: Find income change (I_hat) ---
    # From demand theory, homogeneity implies: ε_YX + ε_YY + η_I_Y = 0
    epsilon_YX = -epsilon_YY - eta_I_Y
    # Change in demand for Y: d_Y_hat = ε_YX * P_X_hat + ε_YY * P_Y_hat + η_I_Y * I_hat
    I_hat = (d_Y_hat - epsilon_YX * P_X_hat - epsilon_YY * P_Y_hat) / eta_I_Y
    print(f"The percentage change in income (I_hat) is derived from the demand for good Y:")
    print(f"d_Y_hat = ε_YX * P_X_hat + ε_YY * P_Y_hat + η_I_Y * I_hat")
    print(f"{d_Y_hat * 100:.2f}% = {epsilon_YX:.2f} * {P_X_hat * 100:.2f}% + {epsilon_YY:.2f} * {P_Y_hat * 100:.2f}% + {eta_I_Y:.2f} * I_hat")
    print(f"Solving for I_hat gives: {I_hat * 100:.2f}%\n")

    # --- Step 5: Find consumption change of good X (d_X_hat) ---
    # From demand theory, homogeneity implies: ε_XX + ε_XY + η_I_X = 0
    epsilon_XY = -epsilon_XX - eta_I_X
    # Change in demand for X: d_X_hat = ε_XX * P_X_hat + ε_XY * P_Y_hat + η_I_X * I_hat
    d_X_hat = epsilon_XX * P_X_hat + epsilon_XY * P_Y_hat + eta_I_X * I_hat
    print(f"The percentage change in consumption of good X (d_X_hat) is derived from its demand function:")
    print(f"d_X_hat = ε_XX * P_X_hat + ε_XY * P_Y_hat + η_I_X * I_hat")
    print(f"d_X_hat = {epsilon_XX:.2f} * {P_X_hat * 100:.2f}% + {epsilon_XY:.2f} * {P_Y_hat * 100:.2f}% + {eta_I_X:.2f} * {I_hat * 100:.2f}%")
    print(f"Solving for d_X_hat gives: {d_X_hat * 100:.2f}%\n")

    # --- Step 6: Format and print the final answer ---
    w_hat_percent = w_hat * 100
    P_X_hat_percent = P_X_hat * 100
    d_X_hat_percent = d_X_hat * 100

    print("Final answer as 3 comma-separated percentage changes (nominal wage, price of X, consumption of X):")
    final_answer_string = f"{w_hat_percent},{P_X_hat_percent},{d_X_hat_percent}"
    print(final_answer_string)
    
    # Return final answer in the specified format
    print(f"<<<{final_answer_string}>>>")


solve_economy_shock()