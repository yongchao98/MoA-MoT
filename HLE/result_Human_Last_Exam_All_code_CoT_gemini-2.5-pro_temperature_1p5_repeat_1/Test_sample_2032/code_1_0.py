import sympy as sp

def solve_variance():
    """
    This function calculates the variance of Y using symbolic integration.
    """
    # Define symbols
    x, u = sp.symbols('x u', real=True)

    # We use symmetry to simplify the problem.
    # We first calculate conditional expectations E[Y|X1=x] and E[Y^2|X1=x].
    # By symmetry of U[0,1], E[Y_x] + E[Y_{1-x}] = 1, which leads to E[Y] = 1/2.
    E_Y = sp.S(1)/2
    print(f"Step 1: Determine E[Y]")
    print(f"By symmetry, E[Y] = {E_Y}\n")

    # The main task is to compute E[Y^2] = integral(E[Y_x^2], (x, 0, 1)).
    # We can show E[Y^2] = integral(1 - 2*E[Y_x] + 2*E[Y_x^2], (x, 0, 1/2)).
    # We need to find expressions for E[Y_x] and E[Y_x^2] for x in [0, 1/2].

    # The formula for the conditional moment is E[Y_x^k] = 6 * integral(u^k * p(u,x) * (1-p(u,x)), (u, 0, 1))
    # For x in [0, 1/2], p(u,x) is a piecewise function:
    # p(u,x) = 2*|u-x| for u in [0, 2x]
    # p(u,x) = u         for u in (2x, 1]

    # Let's define the piecewise components for p(u,x) and the integrand.
    p1 = 2 * (x - u)  # for u in [0, x]
    p2 = 2 * (u - x)  # for u in [x, 2x]
    p3 = u            # for u in (2x, 1]

    # Calculate E[Y_x] for x in [0, 1/2]
    I1_k1 = sp.integrate(u * p1 * (1-p1), (u, 0, x))
    I2_k1 = sp.integrate(u * p2 * (1-p2), (u, x, 2*x))
    I3_k1 = sp.integrate(u * p3 * (1-p3), (u, 2*x, 1))
    E_Yx = sp.simplify(6 * (I1_k1 + I2_k1 + I3_k1))
    
    print(f"Step 2: Calculate the conditional expectation E[Y|X1=x] for x in [0, 1/2]")
    print(f"E[Y_x] = {E_Yx}\n")
    
    # Calculate E[Y_x^2] for x in [0, 1/2]
    I1_k2 = sp.integrate(u**2 * p1 * (1-p1), (u, 0, x))
    I2_k2 = sp.integrate(u**2 * p2 * (1-p2), (u, x, 2*x))
    I3_k2 = sp.integrate(u**2 * p3 * (1-p3), (u, 2*x, 1))
    E_Yx2 = sp.simplify(6 * (I1_k2 + I2_k2 + I3_k2))

    print(f"Step 3: Calculate the conditional second moment E[Y^2|X1=x] for x in [0, 1/2]")
    print(f"E[Y_x^2] = {E_Yx2}\n")

    # Calculate E[Y^2] using the symmetry relation
    integrand_E_Y2 = 1 - 2*E_Yx + 2*E_Yx2
    E_Y2 = sp.integrate(integrand_E_Y2, (x, 0, sp.S(1)/2))
    
    print(f"Step 4: Calculate the second moment E[Y^2]")
    print(f"E[Y^2] = integral(1 - 2*E[Y_x] + 2*E[Y_x^2] dx, (x, 0, 1/2)) = {E_Y2}\n")
    
    # Calculate the variance
    Var_Y = E_Y2 - E_Y**2

    print(f"Step 5: Calculate the variance Var(Y)")
    print(f"Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"Var(Y) = {E_Y2} - ({E_Y})^2")
    print(f"Var(Y) = {E_Y2} - {E_Y**2}")
    print(f"Var(Y) = {Var_Y}")

solve_variance()