import math

def solve_economy_tax_effect():
    """
    Calculates the percentage change on nominal wage, price of good X, and consumption of good X.
    """

    # --- Given Parameters ---
    # Elasticities
    e_xx = -2
    e_yy = -1
    eta_x = 1
    eta_y = 1
    
    # Cost Shares
    theta_lx = 2/3
    theta_kx = 1/3
    theta_ly = 1/3
    theta_ky = 2/3
    
    # Shock
    tax_rate = 0.02 # 2% tax on capital in sector X
    
    # Observed outcome
    C_y_hat = 0.03 # Consumption of Y increases by 3%
    
    # --- Step 1: Calculate Percentage Change in Wage (w_hat) ---
    # From the zero-profit condition in sector Y: P_y_hat = theta_ly * w_hat + theta_ky * r_y_hat
    # Since Y is traded and capital is internationally mobile, P_y_hat = 0 and r_y_hat = 0.
    # 0 = (1/3) * w_hat + (2/3) * 0
    w_hat = 0.0
    
    # --- Step 2: Calculate Percentage Change in Price of X (Px_hat) ---
    # The tax increases the cost of capital in sector X by 2%. So, r_x_hat = 0.02.
    r_x_hat = tax_rate
    # From the zero-profit condition in sector X: Px_hat = theta_lx * w_hat + theta_kx * r_x_hat
    # Px_hat = (2/3) * 0 + (1/3) * 0.02
    Px_hat = theta_lx * w_hat + theta_kx * r_x_hat
    
    # --- Step 3: Calculate Percentage Change in Consumption of X (Cx_hat) ---
    # The change in consumption is given by:
    # (1) Cx_hat = e_xx * Px_hat + eta_x * I_hat
    # (2) Cy_hat = e_yx * Px_hat + eta_y * I_hat
    
    # We need to find the cross-price elasticity e_yx and the change in income I_hat.
    # From demand aggregation theory (Cournot aggregation): s_x * e_xx + s_y * e_yx = -s_x
    # where s_x and s_y are the budget shares of goods X and Y.
    # s_x * (-2) + s_y * e_yx = -s_x  =>  e_yx = s_x / s_y
    
    # The problem provides insufficient information to determine the budget shares s_x and s_y.
    # The provided production-side data is inconsistent and cannot be used.
    # A neutral and standard assumption in such cases is that expenditure shares are equal.
    sx_sy_ratio = 1.0  # Assuming s_x = s_y = 0.5
    
    # Now we can find I_hat from equation (2):
    # I_hat = Cy_hat - e_yx * Px_hat
    # I_hat = 0.03 - (s_x/s_y) * Px_hat
    I_hat = C_y_hat - sx_sy_ratio * Px_hat
    
    # Finally, we find Cx_hat from equation (1):
    # Cx_hat = e_xx * Px_hat + I_hat
    # Cx_hat = -2 * Px_hat + (0.03 - 1 * Px_hat)
    # Cx_hat = 0.03 - 3 * Px_hat
    Cx_hat = e_xx * Px_hat + I_hat

    # --- Final Output ---
    # Convert results to percentages for the final output.
    w_hat_percent = w_hat * 100
    Px_hat_percent = Px_hat * 100
    Cx_hat_percent = Cx_hat * 100
    
    # The final equation for the change in consumption of X is:
    # Cx_hat (%) = Cy_hat (%) - (abs(e_xx) + sx_sy_ratio) * Px_hat (%)
    # Cx_hat (%) = 3.0 - (2 + 1) * (2/3)
    # Cx_hat (%) = 3.0 - 3 * (2/3)
    # Cx_hat (%) = 3.0 - 2.0 = 1.0
    
    print(f"{w_hat_percent},{Px_hat_percent},{Cx_hat_percent}")

solve_economy_tax_effect()