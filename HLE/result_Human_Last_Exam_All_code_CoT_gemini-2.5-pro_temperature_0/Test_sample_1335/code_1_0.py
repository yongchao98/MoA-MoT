import sys

# This is a helper function to format the output as requested.
def print_equation_step(step_number, description, equation_str, values, result_name, result_value, unit="%"):
    """Helper function to print calculation steps."""
    # Replace placeholders in the equation string with formatted values
    formatted_equation = equation_str
    for i, val in enumerate(values):
        # Format floats with precision, keep integers as is
        if isinstance(val, float):
            formatted_val = f"{val:.4f}"
        else:
            formatted_val = str(val)
        formatted_equation = formatted_equation.replace(f"v{i+1}", formatted_val)

    print(f"Step {step_number}: {description}")
    print(f"Equation: {formatted_equation}")
    print(f"Result: {result_name} = {result_value * 100:.4f}{unit}")
    print("-" * 30)

def solve_economy_problem():
    """
    Solves the economic problem based on the provided parameters.
    """
    # --- Given Parameters ---
    tax_on_capital_x = 0.02  # 2% tax
    consumption_y_change = 0.03 # 3% increase

    # Cost shares
    theta_LX = 2/3
    theta_KX = 1/3
    theta_LY = 1/3
    
    # Demand elasticities
    eta_XX = -2.0
    eta_YY = -1.0
    epsilon_XI = 1.0
    epsilon_YI = 1.0

    # --- Step 1: Calculate Percentage Change in Wage (w_hat) ---
    # From zero-profit condition in sector Y: P_hat_Y = theta_LY * w_hat + theta_KY * r_hat_Y
    # P_hat_Y = 0 (traded good), r_hat_Y = 0 (world rate, no tax in Y)
    # 0 = theta_LY * w_hat + 0 => w_hat = 0
    w_hat = 0.0
    print_equation_step(
        step_number=1,
        description="Calculate percentage change in nominal wage (w_hat)",
        equation_str="0 = v1 * w_hat",
        values=[theta_LY],
        result_name="w_hat",
        result_value=w_hat
    )

    # --- Step 2: Calculate Percentage Change in Price of X (Px_hat) ---
    # From zero-profit condition in sector X: Px_hat = theta_LX * w_hat + theta_KX * r_hat_X
    # r_hat_X is the change in the cost of capital due to the tax
    r_hat_X = tax_on_capital_x
    Px_hat = theta_LX * w_hat + theta_KX * r_hat_X
    print_equation_step(
        step_number=2,
        description="Calculate percentage change in price of good X (Px_hat)",
        equation_str="Px_hat = v1 * v2 + v3 * v4",
        values=[theta_LX, w_hat, theta_KX, r_hat_X],
        result_name="Px_hat",
        result_value=Px_hat
    )

    # --- Step 3: Calculate Percentage Change in Income (I_hat) ---
    # Use the given information about consumption of Y.
    # C_hat_Y = eta_YX * Px_hat + eta_YY * Py_hat + epsilon_YI * I_hat
    # First, find eta_YX from homogeneity: eta_YX + eta_YY + epsilon_YI = 0
    eta_YX = -eta_YY - epsilon_YI
    # Py_hat = 0
    # I_hat = (C_hat_Y - eta_YX * Px_hat) / epsilon_YI
    I_hat = (consumption_y_change - eta_YX * Px_hat) / epsilon_YI
    print("Step 3: Calculate percentage change in income (I_hat)")
    print(f"From demand homogeneity for Y: eta_YX + eta_YY + epsilon_YI = 0 => eta_YX = {-eta_YY - epsilon_YI}")
    print_equation_step(
        step_number=3.1,
        description="Using the demand equation for Y",
        equation_str="v1 = v2 * v3 + v4 * 0 + v5 * I_hat",
        values=[consumption_y_change, eta_YX, Px_hat, eta_YY, epsilon_YI],
        result_name="I_hat",
        result_value=I_hat
    )

    # --- Step 4: Calculate Percentage Change in Consumption of X (Cx_hat) ---
    # C_hat_X = eta_XX * Px_hat + eta_XY * Py_hat + epsilon_XI * I_hat
    # First, find eta_XY from homogeneity: eta_XX + eta_XY + epsilon_XI = 0
    eta_XY = -eta_XX - epsilon_XI
    # Py_hat = 0
    Cx_hat = eta_XX * Px_hat + eta_XY * 0 + epsilon_XI * I_hat
    print("Step 4: Calculate percentage change in consumption of X (Cx_hat)")
    print(f"From demand homogeneity for X: eta_XX + eta_XY + epsilon_XI = 0 => eta_XY = {-eta_XX - epsilon_XI}")
    print_equation_step(
        step_number=4.1,
        description="Using the demand equation for X",
        equation_str="Cx_hat = v1 * v2 + v3 * 0 + v4 * v5",
        values=[eta_XX, Px_hat, eta_XY, epsilon_XI, I_hat],
        result_name="Cx_hat",
        result_value=Cx_hat
    )

    # --- Final Answer ---
    # The question asks for percentage change on nominal wage, price of good X, and the consumption of good X
    final_w_hat_percent = w_hat * 100
    final_Px_hat_percent = Px_hat * 100
    final_Cx_hat_percent = Cx_hat * 100

    # Print the final answer in the required format
    print(f"{final_w_hat_percent},{final_Px_hat_percent},{final_Cx_hat_percent}")

solve_economy_problem()