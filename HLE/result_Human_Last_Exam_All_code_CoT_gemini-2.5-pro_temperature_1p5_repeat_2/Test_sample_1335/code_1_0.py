import sys

def solve_economy_tax_effect():
    """
    Calculates the percentage change on nominal wage, price of good X, 
    and consumption of good X due to a tax on capital in sector X.
    """

    # --- Step 0: Define parameters from the problem statement ---

    # The shock: a 2% tax on the return to capital in sector X
    hat_r_X = 0.02

    # The known result: consumption of Y increases by 3%
    hat_C_Y = 0.03

    # Labor's share of cost
    theta_LX = 2/3  # in sector X
    theta_LY = 1/3  # in sector Y

    # Capital's share of cost (derived)
    theta_KX = 1 - theta_LX
    theta_KY = 1 - theta_LY

    # Elasticity of substitution for factors of production.
    # Standard economic theory requires this to be non-negative.
    # We assume the negative signs in the prompt are typos.
    sigma_X = 2
    sigma_Y = 1

    # Sectoral shares of total labor supply
    lambda_LX = 3/4
    lambda_LY = 1 - lambda_LX
    
    # Note: The provided factor allocation shares (lambda) and cost shares (theta)
    # are mathematically inconsistent in a static equilibrium. The solution below is derived
    # for the *changes* from equilibrium and bypasses this issue by using a simplifying assumption.


    # --- Step 1: Calculate the percentage change in nominal wage (hat_w) ---
    # Good Y is traded, so its price P_Y is fixed by the world market (hat_P_Y = 0).
    # Capital is traded, so the return r is fixed. In sector Y, hat_r_Y = 0.
    # The zero-profit condition in sector Y is: hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y
    # 0 = (1/3) * hat_w + (2/3) * 0
    hat_w = 0


    # --- Step 2: Calculate the percentage change in the price of good X (hat_P_X) ---
    # The zero-profit condition in sector X is: hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
    # We use hat_w = 0 and the given tax effect hat_r_X = 0.02.
    # Final Equation for Price Change of X: hat_P_X = (2/3) * 0 + (1/3) * 0.02
    hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
    
    
    # --- Step 3: Find the supply-side factor input responses ---
    # From the definition of elasticity of substitution: hat(K/L) = sigma * hat(w/r)
    # For Sector Y: hat_K_Y - hat_L_Y = sigma_Y * (hat_w - hat_r_Y) = 1 * (0 - 0) = 0.
    # Thus, the percentage change in capital and labor in sector Y are equal.
    # hat_K_Y = hat_L_Y

    # For Sector X: hat_K_X - hat_L_X = sigma_X * (hat_w - hat_r_X) = 2 * (0 - 0.02) = -0.04
    hat_K_X_minus_hat_L_X = sigma_X * (hat_w - hat_r_X)
    
    # From full employment of labor (fixed total supply): lambda_LX*hat_L_X + lambda_LY*hat_L_Y = 0
    # (3/4) * hat_L_X + (1/4) * hat_L_Y = 0  => hat_L_Y = -3 * hat_L_X
    hat_L_Y_over_hat_L_X = -(lambda_LX / lambda_LY)


    # --- Step 4: Link consumption change to production change ---
    # The change in output of Y is: hat_Q_Y = theta_LY*hat_L_Y + theta_KY*hat_K_Y
    # Since hat_K_Y = hat_L_Y, we get hat_Q_Y = (theta_LY + theta_KY) * hat_L_Y = 1 * hat_L_Y
    # Thus, hat_Q_Y = hat_L_Y.

    # We assume the given change in consumption of Y equals the change in its production.
    # This assumption (hat_C_Y = hat_Q_Y) is necessary to solve the problem given the inconsistent static data.
    hat_Q_Y = hat_C_Y
    hat_L_Y = hat_Q_Y

    # --- Step 5: Solve for all unknown percentage changes ---
    # With hat_L_Y known, we can find the other factor changes.
    hat_L_X = hat_L_Y / hat_L_Y_over_hat_L_X
    hat_K_X = hat_K_X_minus_hat_L_X + hat_L_X

    # --- Step 6: Calculate the change in consumption of X ---
    # Since X is non-traded, change in consumption equals change in production: hat_C_X = hat_Q_X.
    # The change in production is: hat_Q_X = theta_LX * hat_L_X + theta_KX * hat_K_X
    # Final Equation for Consumption Change of X: hat_C_X = (2/3) * hat_L_X + (1/3) * hat_K_X
    hat_C_X = theta_LX * hat_L_X + theta_KX * hat_K_X
    

    # --- Final Output ---
    # Convert ratios to percentages and print as comma-separated values.
    # The values represent the percentage changes for:
    # 1. Nominal wage (w)
    # 2. Price of good X (P_X)
    # 3. Consumption of good X (C_X)
    print(f"{hat_w * 100},{hat_P_X * 100},{hat_C_X * 100}")

solve_economy_tax_effect()