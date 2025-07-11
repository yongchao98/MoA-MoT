from sympy import symbols, solve, Eq

def solve_economy_tax_effect():
    """
    Solves for the percentage changes in wage, price of X, and consumption of X
    based on the given economic model and a tax on capital in sector X.
    """
    
    # --- 1. Define Parameters (in decimal form) ---
    
    # Elasticities
    sigma_X = 2.0  # Elasticity of substitution in X (assuming positive value)
    sigma_Y = 1.0  # Elasticity of substitution in Y (assuming positive value, Cobb-Douglas)
    eta_IX = 1.0   # Income elasticity of demand for X
    eta_IY = 1.0   # Income elasticity of demand for Y
    epsilon_XX = -2.0 # Own-price elasticity of demand for X

    # Factor allocation shares (lambda)
    lambda_LX = 3/4
    lambda_LY = 1/4
    lambda_KX = 1/4
    lambda_KY = 3/4
    
    # Cost shares (theta)
    theta_LX = 2/3
    theta_KX = 1 - theta_LX
    theta_LY = 1/3
    theta_KY = 1 - theta_LY
    
    # The shock: 2% tax on capital return in X
    tax = 0.02
    
    # --- 2. Solve for Factor and Good Price Changes (% change) ---
    
    # External conditions for a small open economy
    # Percentage change in Price of Y is 0 (traded good)
    hat_P_Y = 0
    # Net return to capital is fixed by world rate, so change is 0.
    hat_r_K = 0
    # Producers in Y pay the world rate.
    hat_r_Y = hat_r_K
    # Producers in X pay the world rate plus the tax.
    hat_r_X = hat_r_K + tax
    
    # From price-cost equation for Sector Y: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y
    # 0 = (1/3) * hat_w + (2/3) * 0
    hat_w = 0.0
    
    # From price-cost equation for Sector X: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
    hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
    
    # --- 3. Determine Changes in Sectoral Output ---
    
    # Define symbolic variables for the output changes
    hat_Q_X, hat_Q_Y = symbols('hat_Q_X hat_Q_Y')

    # Change in factor demand equations:
    # dL_i/L_i = dQ_i/Q_i - theta_Ki * sigma_i * (dw/w - dr_i/r_i)
    # dK_i/K_i = dQ_i/Q_i + theta_Li * sigma_i * (dw/w - dr_i/r_i)
    hat_L_X = hat_Q_X - theta_KX * sigma_X * (hat_w - hat_r_X)
    hat_K_X = hat_Q_X + theta_LX * sigma_X * (hat_w - hat_r_X)
    
    hat_L_Y = hat_Q_Y - theta_KY * sigma_Y * (hat_w - hat_r_Y)
    hat_K_Y = hat_Q_Y + theta_LY * sigma_Y * (hat_w - hat_r_Y)
    
    # Full employment equations (fixed total factor supplies)
    # lambda_LX * hat_L_X + lambda_LY * hat_L_Y = 0
    # lambda_KX * hat_K_X + lambda_KY * hat_K_Y = 0
    eq1 = Eq(lambda_LX * hat_L_X + lambda_LY * hat_L_Y, 0)
    eq2 = Eq(lambda_KX * hat_K_X + lambda_KY * hat_K_Y, 0)
    
    # Solve the system of two linear equations for hat_Q_X and hat_Q_Y
    output_changes = solve((eq1, eq2), (hat_Q_X, hat_Q_Y))
    hat_Q_X_val = output_changes[hat_Q_X]

    # --- 4. Solve for Consumption of Good X ---
    
    # For a non-traded good, change in consumption = change in production
    hat_D_X = hat_Q_X_val

    # Convert results to percentages for output
    wage_change_pct = hat_w * 100
    price_X_change_pct = hat_P_X * 100
    consumption_X_change_pct = hat_D_X * 100
    
    # Print the results as comma-separated values
    print("Percentage change in nominal wage:", wage_change_pct)
    print("Percentage change in price of good X:", price_X_change_pct)
    print("Percentage change in consumption of good X:", consumption_X_change_pct)
    
    # Final answer in requested format
    final_answer_string = f"{wage_change_pct},{price_X_change_pct},{consumption_X_change_pct}"
    return final_answer_string

# Execute the function and print the final answer
final_answer = solve_economy_tax_effect()
print(f"<<<{final_answer}>>>")