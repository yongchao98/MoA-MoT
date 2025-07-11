import sys
from io import StringIO

def solve_economy_tax_effect():
    """
    This function calculates the percentage changes in nominal wage, price of good X,
    and consumption of good X based on the given economic model and a tax shock.
    """
    # --- 1. Define known parameters and shocks ---
    # Cost shares
    theta_LX = 2/3  # Labor's share of cost in sector X
    theta_KX = 1/3  # Capital's share of cost in sector X
    theta_LY = 1/3  # Labor's share of cost in sector Y
    theta_KY = 2/3  # Capital's share of cost in sector Y

    # Demand elasticities
    epsilon_XX = -2.0  # Own price elasticity of demand for good X
    epsilon_YY = -1.0  # Own price elasticity of demand for good Y
    eta_IX = 1.0     # Income elasticity of demand for both goods
    eta_IY = 1.0

    # Initial shocks and fixed values
    r_X_hat = 0.02   # Percentage change in cost of capital for sector X (due to 2% tax)
    P_Y_hat = 0.0    # Good Y is freely traded, so its price is fixed
    r_Y_hat = 0.0    # Capital cost for sector Y is the unchanged world rate
    D_Y_hat = 0.03   # Given consumption of Y increases by 3%

    # --- 2. Calculate Percentage Change in Nominal Wage (w_hat) ---
    # From the zero-profit condition in sector Y: P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat
    # 0 = (1/3) * w_hat + (2/3) * 0  => w_hat = 0
    w_hat = (P_Y_hat - theta_KY * r_Y_hat) / theta_LY
    
    print("Step 1: Calculate percentage change in nominal wage (w_hat)")
    print(f"Equation: P_Y_hat = theta_LY * w_hat + theta_KY * r_Y_hat")
    print(f"{P_Y_hat:.2f} = {theta_LY:.2f} * w_hat + {theta_KY:.2f} * {r_Y_hat:.2f}")
    print(f"Result: w_hat = {w_hat:.4f}\n")

    # --- 3. Calculate Percentage Change in Price of Good X (P_X_hat) ---
    # From the zero-profit condition in sector X: P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
    P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat
    
    print("Step 2: Calculate percentage change in price of good X (P_X_hat)")
    print(f"Equation: P_X_hat = theta_LX * w_hat + theta_KX * r_X_hat")
    print(f"P_X_hat = {theta_LX:.2f} * {w_hat:.2f} + {theta_KX:.2f} * {r_X_hat:.2f}")
    print(f"Result: P_X_hat = {P_X_hat:.4f}\n")
    
    # --- 4. Calculate Percentage Change in Nominal Income (I_hat) ---
    # From demand function for Y: D_Y_hat = epsilon_YX * P_X_hat + epsilon_YY * P_Y_hat + eta_IY * I_hat
    # First find epsilon_YX from homogeneity: epsilon_YX + epsilon_YY + eta_IY = 0
    epsilon_YX = - epsilon_YY - eta_IY
    # Now solve for I_hat
    I_hat = (D_Y_hat - epsilon_YX * P_X_hat - epsilon_YY * P_Y_hat) / eta_IY

    print("Step 3: Calculate percentage change in nominal income (I_hat)")
    print(f"From demand homogeneity for Y: epsilon_YX + epsilon_YY + eta_IY = 0")
    print(f"epsilon_YX + ({epsilon_YY:.2f}) + {eta_IY:.2f} = 0  => epsilon_YX = {-epsilon_YY - eta_IY:.2f}")
    print(f"From demand function for Y: D_Y_hat = epsilon_YX * P_X_hat + epsilon_YY * P_Y_hat + eta_IY * I_hat")
    print(f"{D_Y_hat:.2f} = {epsilon_YX:.2f} * {P_X_hat:.4f} + {epsilon_YY:.2f} * {P_Y_hat:.2f} + {eta_IY:.2f} * I_hat")
    print(f"Result: I_hat = {I_hat:.4f}\n")

    # --- 5. Calculate Percentage Change in Consumption of Good X (D_X_hat) ---
    # From demand function for X: D_X_hat = epsilon_XX * P_X_hat + epsilon_XY * P_Y_hat + eta_IX * I_hat
    # With P_Y_hat = 0, this simplifies to D_X_hat = epsilon_XX * P_X_hat + eta_IX * I_hat
    D_X_hat = epsilon_XX * P_X_hat + eta_IX * I_hat

    print("Step 4: Calculate percentage change in consumption of good X (D_X_hat)")
    print(f"Equation: D_X_hat = epsilon_XX * P_X_hat + eta_IX * I_hat")
    print(f"D_X_hat = {epsilon_XX:.2f} * {P_X_hat:.4f} + {eta_IX:.2f} * {I_hat:.4f}")
    print(f"Result: D_X_hat = {D_X_hat:.4f}\n")
    
    # --- 6. Final Answer ---
    # The question asks for the percentage changes
    w_hat_percent = w_hat * 100
    P_X_hat_percent = P_X_hat * 100
    D_X_hat_percent = D_X_hat * 100
    
    print("Final percentage changes (wage, price of X, consumption of X):")
    # Using .15f to show the repeating decimals accurately
    print(f"{w_hat_percent:.15f},{P_X_hat_percent:.15f},{D_X_hat_percent:.15f}")
    
    return w_hat_percent, P_X_hat_percent, D_X_hat_percent

# Execute the function and capture output
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

w_p, p_p, d_p = solve_economy_tax_effect()

sys.stdout = old_stdout
print(captured_output.getvalue())

# The final answer in the required format
final_answer_str = f"{w_p},{p_p},{d_p}"
print(f"<<<{final_answer_str}>>>")