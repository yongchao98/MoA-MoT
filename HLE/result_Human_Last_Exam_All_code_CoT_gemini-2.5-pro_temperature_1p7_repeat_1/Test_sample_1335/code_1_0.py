def solve_economy_shock():
    """
    Calculates the percentage changes in wage, price of X, and consumption of X.
    """
    # Parameters from the problem statement (in percentage terms)
    tax_on_capital_x = 2.0
    consumption_y_change = 3.0

    # Factor cost shares
    theta_Lx = 2/3
    theta_Kx = 1 - theta_Lx
    theta_Ly = 1/3
    theta_Ky = 1 - theta_Ly

    # Demand elasticities
    own_price_elast_x = -2.0
    own_price_elast_y = -1.0
    income_elast_x = 1.0
    income_elast_y = 1.0

    # ----- Step 1: Calculate the percentage change in the nominal wage (w_hat) -----
    # In sector Y, P_y_hat = theta_Ly * w_hat + theta_Ky * r_y_hat.
    # Since Y is traded, P_y is fixed, so P_y_hat = 0.
    # Since capital is mobile, the return r in sector Y is fixed, so r_y_hat = 0.
    # The equation is: 0 = (1/3) * w_hat + (2/3) * 0, which implies w_hat = 0.
    w_hat = 0.0

    # ----- Step 2: Calculate the percentage change in the price of good X (Px_hat) -----
    # In sector X, Px_hat = theta_Lx * w_hat + theta_Kx * r_x_hat.
    # The 2% tax means r_x_hat = 2.0. We found w_hat = 0.
    # The equation is: Px_hat = (2/3) * 0 + (1/3) * 2.0.
    r_x_hat = tax_on_capital_x
    Px_hat = theta_Lx * w_hat + theta_Kx * r_x_hat

    # ----- Step 3: Calculate the percentage change in the consumption of good X (Cx_hat) -----
    # First, find the change in income (I_hat). From demand theory for good Y:
    # Cy_hat = elasticity_yx * Px_hat + elasticity_yy * Py_hat + income_elasticity_y * I_hat.
    # Homogeneity (e_yx + e_yy + n_y = 0) implies e_yx + (-1) + 1 = 0, so e_yx = 0.
    # With Py_hat = 0, the equation simplifies to Cy_hat = income_elasticity_y * I_hat.
    # Since income_elasticity_y = 1, Cy_hat = I_hat.
    # Given Cy_hat = 3.0, it follows that I_hat = 3.0.
    I_hat = consumption_y_change

    # Now use the demand function for X:
    # Cx_hat = own_price_elast_x * Px_hat + income_elast_x * I_hat
    # The equation is: Cx_hat = -2.0 * Px_hat + 1.0 * I_hat.
    Cx_hat = own_price_elast_x * Px_hat + income_elast_x * I_hat

    # ----- Print results with equations -----

    print("The percentage changes are calculated as follows (values are in percent):\n")

    print(f"1. Change in nominal wage (w_hat):")
    print(f"   Equation: 0 = ({theta_Ly:.3f}) * w_hat + ({theta_Ky:.3f}) * 0")
    print(f"   w_hat = {w_hat:.4f}%\n")

    print(f"2. Change in price of good X (Px_hat):")
    print(f"   Equation: Px_hat = ({theta_Lx:.3f}) * ({w_hat:.2f}) + ({theta_Kx:.3f}) * ({r_x_hat:.2f})")
    print(f"   Px_hat = {Px_hat:.4f}%\n")

    print(f"3. Change in consumption of good X (Cx_hat):")
    print(f"   First, income change I_hat = Cy_hat = {I_hat:.2f}%")
    print(f"   Equation: Cx_hat = ({own_price_elast_x:.2f}) * ({Px_hat:.4f}) + ({income_elast_x:.2f}) * ({I_hat:.2f})")
    print(f"   Cx_hat = {Cx_hat:.4f}%\n")

    print("Final answers as 3 comma-separated values (percentage change in nominal wage, price of good X, consumption of good X):")
    print(f"{w_hat},{Px_hat},{Cx_hat}")

solve_economy_shock()