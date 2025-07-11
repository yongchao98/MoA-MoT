import sys

def solve_economy_tax_effect():
    """
    Calculates the percentage change on nominal wage, price of good X,
    and consumption of good X due to a tax on capital in sector X.
    """

    # 1. Define the given parameters from the problem description.
    # Note: "hat" variables represent proportional changes (e.g., 0.02 for 2%).

    # Tax on return to capital in sector X
    hat_r_X = 0.02

    # Labor's share of cost in sector X and Y
    theta_LX = 2/3
    theta_LY = 1/3

    # Capital's share of cost in sector X and Y
    theta_KX = 1/3
    theta_KY = 2/3

    # Change in price of good Y (traded good, price is fixed)
    hat_P_Y = 0

    # Change in cost of capital for sector Y (world rate is fixed)
    hat_r_Y = 0

    # Given change in consumption of Y
    hat_D_Y = 0.03

    # Demand elasticities for good Y
    eps_YY = -1.0  # Own-price elasticity
    eta_Y = 1.0    # Income elasticity

    # Demand elasticities for good X
    eps_XX = -2.0  # Own-price elasticity
    eta_X = 1.0    # Income elasticity

    # 2. Calculate the percentage change in nominal wage (hat_w).
    # From the zero-profit condition in the traded goods sector (Y):
    # hat_P_Y = theta_LY * hat_w + theta_KY * hat_r_Y
    # 0 = (1/3) * hat_w + (2/3) * 0
    # This implies hat_w = 0.
    hat_w = 0.0

    # 3. Calculate the percentage change in the price of good X (hat_P_X).
    # From the zero-profit condition in sector X:
    # hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X
    hat_P_X = theta_LX * hat_w + theta_KX * hat_r_X

    # 4. Calculate the percentage change in nominal income (hat_I).
    # First, find the cross-price elasticity of demand for Y w.r.t. price of X (eps_YX).
    # From the homogeneity of degree zero property of demand functions: eps_YX + eps_YY + eta_Y = 0
    eps_YX = -eps_YY - eta_Y

    # Next, use the demand equation for Y to solve for hat_I:
    # hat_D_Y = eps_YX * hat_P_X + eps_YY * hat_P_Y + eta_Y * hat_I
    # We can rearrange this to solve for hat_I:
    # hat_I = (hat_D_Y - eps_YX * hat_P_X - eps_YY * hat_P_Y) / eta_Y
    hat_I = (hat_D_Y - eps_YX * hat_P_X - eps_YY * hat_P_Y) / eta_Y

    # 5. Calculate the percentage change in consumption of good X (hat_D_X).
    # Use the demand equation for X:
    # hat_D_X = eps_XX * hat_P_X + (cross-price term) + eta_X * hat_I
    # The cross-price term involves hat_P_Y, which is 0, so the term is 0.
    hat_D_X = eps_XX * hat_P_X + eta_X * hat_I

    # 6. Convert results to percentages for the final output.
    wage_change_perc = hat_w * 100
    price_X_change_perc = hat_P_X * 100
    consumption_X_change_perc = hat_D_X * 100

    # Print the final answer as 3 comma-separated values.
    print(f"{wage_change_perc},{price_X_change_perc},{consumption_X_change_perc}")

solve_economy_tax_effect()
<<<0.0,0.6666666666666666,1.6666666666666667>>>