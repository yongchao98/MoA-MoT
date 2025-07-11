import math

def solve_economy_tax_effect():
    """
    This function calculates the percentage changes in nominal wage, the price of good X,
    and the consumption of good X based on the given economic model and a tax shock.
    """
    
    # 1. Define given parameters and shocks
    
    # Shock: 2% tax on return to capital in sector X
    # This increases the cost of capital for firms in sector X by 2%
    hat_rX = 0.02
    
    # Factor cost shares
    theta_LX = 2/3  # Labor's share of cost in sector X
    theta_LY = 1/3  # Labor's share of cost in sector Y
    theta_KX = 1 - theta_LX  # Capital's share of cost in sector X
    theta_KY = 1 - theta_LY  # Capital's share of cost in sector Y
    
    # Demand parameters
    eta_X = 1.0  # Income elasticity of demand for X
    eta_Y = 1.0  # Income elasticity of demand for Y
    epsilon_XX = -2.0  # Own price elasticity of demand for X
    epsilon_YY = -1.0  # Own price elasticity of demand for Y
    
    # Given outcome from the shock
    hat_DY = 0.03  # Consumption of Y increases by 3%

    # Model assumptions for small open economy
    # Price of traded good Y is fixed, so its percentage change is 0.
    hat_PY = 0.0
    # Cost of capital for firms in Y is the world rental rate, which is fixed.
    hat_rY = 0.0
    
    # 2. Calculate the percentage change in nominal wage (hat_w)
    # The zero-profit condition for sector Y is: hat_PY = theta_LY * hat_w + theta_KY * hat_rY
    # Rearranging for hat_w: hat_w = (hat_PY - theta_KY * hat_rY) / theta_LY
    # Substituting values: hat_w = (0 - (2/3) * 0) / (1/3) = 0
    hat_w = (hat_PY - theta_KY * hat_rY) / theta_LY
    
    # 3. Calculate the percentage change in the price of good X (hat_PX)
    # The zero-profit condition for sector X is: hat_PX = theta_LX * hat_w + theta_KX * hat_rX
    hat_PX = theta_LX * hat_w + theta_KX * hat_rX
    
    # 4. Calculate the percentage change in the consumption of good X (hat_DX)
    
    # First, find the percentage change in consumer income (hat_I_cons)
    # The demand equation for Y is: hat_DY = epsilon_YX * hat_PX + epsilon_YY * hat_PY + eta_Y * hat_I_cons
    # From homogeneity of demand for Y: epsilon_YX + epsilon_YY + eta_Y = 0
    epsilon_YX = -epsilon_YY - eta_Y # This gives 0
    # With hat_PY = 0, the equation simplifies to: hat_DY = epsilon_YX * hat_PX + eta_Y * hat_I_cons
    # Rearranging for hat_I_cons: hat_I_cons = (hat_DY - epsilon_YX * hat_PX) / eta_Y
    hat_I_cons = (hat_DY - epsilon_YX * hat_PX) / eta_Y

    # Next, find the percentage change in consumption of good X
    # The demand equation for X is: hat_DX = epsilon_XX * hat_PX + epsilon_XY * hat_PY + eta_X * hat_I_cons
    # We only need to plug in the values as hat_PY is 0
    hat_DX = epsilon_XX * hat_PX + eta_X * hat_I_cons
    
    # 5. Format and print results
    
    # Convert decimal changes to percentages
    pct_w = hat_w * 100
    pct_PX = hat_PX * 100
    pct_DX = hat_DX * 100

    print("---Calculations---")
    print(f"Percentage change in wage (w):")
    print(f"Formula: (%ΔPY - θ_KY * %Δr_Y) / θ_LY")
    print(f"({hat_PY*100} - {round(theta_KY, 3)} * {hat_rY*100}) / {round(theta_LY, 3)} = {pct_w}")

    print(f"\nPercentage change in price of X (PX):")
    print(f"Formula: θ_LX * %Δw + θ_KX * %ΔrX")
    print(f"{round(theta_LX, 3)} * {pct_w} + {round(theta_KX, 3)} * {hat_rX*100} = {pct_PX}")

    print(f"\nPercentage change in consumption of X (DX):")
    # To find this, we first needed %ΔI = (%ΔDY - ε_YX * %ΔPX) / η_Y
    # %ΔI = (3 - 0 * {pct_PX}) / 1.0 = 3.0
    pct_I_cons = hat_I_cons * 100
    print(f"Formula: ε_XX * %ΔPX + η_X * %ΔI_cons")
    print(f"{epsilon_XX} * {pct_PX} + {eta_X} * {pct_I_cons} = {pct_DX}")
    print("---End of Calculations---")
    
    # Print the final answer in the required comma-separated format
    final_answer = f"{pct_w},{pct_PX},{pct_DX}"
    print(f"\nFinal comma-separated values:")
    print(final_answer)
    
    # Return final answer in specified format
    return f"<<<{pct_w},{pct_PX},{pct_DX}>>>"

# Execute the function and capture the final answer for the 'answer content' format.
final_answer_formatted = solve_economy_tax_effect()
# This final print is for the required format <<<answer>>> at the end of the entire response.
# print(final_answer_formatted)