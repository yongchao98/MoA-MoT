import math

def solve_economy_tax_effect():
    """
    This function calculates the percentage change in nominal wage, price of good X,
    and consumption of good X based on the given economic model and a tax shock.
    """
    
    # --- Step 0: Define constants from the problem statement ---
    # Note: All "_change" variables represent percentage changes in decimal form (e.g., 2% = 0.02)
    
    # Tax on return to capital in sector X
    tax_rate_change = 0.02
    
    # Given resulting change in consumption of Y
    consumption_Y_change = 0.03

    # Income elasticities of demand
    eta_IX = 1.0  # for X
    eta_IY = 1.0  # for Y

    # Own-price elasticities of demand
    eta_XX = -2.0 # for X
    eta_YY = -1.0 # for Y

    # Factor cost shares
    theta_LX = 2.0 / 3.0 # Labor's share in X
    theta_KX = 1.0 / 3.0 # Capital's share in X
    theta_LY = 1.0 / 3.0 # Labor's share in Y
    theta_KY = 2.0 / 3.0 # Capital's share in Y

    # Changes for traded goods and factors
    # Price of Y is fixed by world market
    P_Y_change = 0.0
    # Cost of capital in sector Y is the world rate r, which is constant
    r_Y_change = 0.0
    # Cost of capital in sector X is r*(1+t), so its percentage change is the tax itself
    r_X_change = tax_rate_change

    # --- Step 1: Calculate percentage change in nominal wage (w_change) ---
    # From the zero-profit condition in the traded goods sector Y:
    # P_Y_change = theta_LY * w_change + theta_KY * r_Y_change
    # 0 = (1/3) * w_change + (2/3) * 0
    # This implies w_change = 0
    w_change = 0.0
    
    # --- Step 2: Calculate percentage change in the price of good X (P_X_change) ---
    # From the zero-profit condition in sector X:
    # P_X_change = theta_LX * w_change + theta_KX * r_X_change
    P_X_change = theta_LX * w_change + theta_KX * r_X_change
    
    # --- Step 3: Calculate percentage change in the consumption of good X (C_X_change) ---
    # Part A: Find the change in nominal income (I_change) using data from good Y
    # The demand function change for Y is:
    # C_Y_change = eta_YX * P_X_change + eta_YY * P_Y_change + eta_IY * I_change
    # We find the cross-price elasticity eta_YX from the homogeneity condition:
    # eta_YX + eta_YY + eta_IY = 0  => eta_YX - 1.0 + 1.0 = 0  => eta_YX = 0
    eta_YX = 0.0
    # Now solve for I_change:
    # 0.03 = (0 * P_X_change) + (-1.0 * 0) + (1.0 * I_change)
    I_change = consumption_Y_change / eta_IY
    
    # Part B: Find the change in consumption of X
    # The demand function change for X is:
    # C_X_change = eta_XX * P_X_change + eta_XY * P_Y_change + eta_IX * I_change
    # We find the cross-price elasticity eta_XY from the homogeneity condition:
    # eta_XX + eta_XY + eta_IX = 0  => -2.0 + eta_XY + 1.0 = 0 => eta_XY = 1.0
    eta_XY = 1.0
    # Now solve for C_X_change:
    C_X_change = eta_XX * P_X_change + eta_XY * P_Y_change + eta_IX * I_change
    
    # --- Final Step: Format and print the results ---
    # The required answers are the percentage changes (not decimals).
    w_change_pct = w_change * 100
    P_X_change_pct = P_X_change * 100
    C_X_change_pct = C_X_change * 100
    
    # Print the final result as 3 comma-separated values.
    # The values correspond to: % change in nominal wage, % change in price of good X, % change in consumption of good X.
    print(f"{w_change_pct},{P_X_change_pct},{C_X_change_pct}")

solve_economy_tax_effect()