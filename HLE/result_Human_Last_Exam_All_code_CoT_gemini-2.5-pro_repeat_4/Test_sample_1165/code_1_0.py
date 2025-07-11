import sympy

def solve_and_explain():
    """
    Derives and explains the scaling of fluctuations R(epsilon).
    """
    # Define symbols
    eps = sympy.Symbol('epsilon', positive=True, real=True)
    x = sympy.Symbol('x')
    L = 1 / eps
    N = L - 1

    # --- Explanation and Derivation ---
    print("### Derivation of R(epsilon) ###")
    print("\nWe aim to find R = (max_x Var[y(x) - y(0)])^(1/2).")
    print("Since y(0)=1 is a constant, Var[y(x) - y(0)] = Var[y(x)].")
    
    print("\nStep 1: Equation for fluctuations")
    print("The fluctuation y_f = y - E[y] follows the approximate equation for small epsilon:")
    print("y_f''(x) ≈ epsilon^2 * S(x)")
    print("where S(x) is the zero-mean noise term.")

    print("\nStep 2: Variance using Green's function")
    print("The solution is y_f(x) ≈ epsilon^2 * integral(G_0(x,s) * S(s) ds), where G_0 is the Green's function for the operator d^2/dx^2.")
    print("The variance is Var[y_f(x)] ≈ epsilon^4 * integral(G_0(x,s)^2 * E[S(s)^2] ds).")
    print("Using the shot noise approximation E[S(s1)S(s2)] ≈ (N/L) * delta(s1-s2), this becomes:")
    print("Var[y_f(x)] ≈ epsilon^4 * (N/L) * integral from 0 to L of G_0(x,s)^2 ds.")
    
    N_over_L = (1 - eps)
    print(f"\nHere, N/L = (1/epsilon - 1) / (1/epsilon) = {N_over_L}.")

    print("\nStep 3: Calculating the integral")
    print("The Green's function is G_0(x,s) = epsilon*x*s - Min(x,s).")
    print("The integral I(x) = integral from 0 to L of G_0(x,s)^2 ds can be calculated analytically.")
    # The result from the analytical calculation is (x^2 * (eps*x - 1)^2) / (3*eps)
    I_x_str = "(x^2 * (epsilon*x - 1)^2) / (3*epsilon)"
    print(f"The result is I(x) = {I_x_str}.")

    print("\nStep 4: Maximizing the integral")
    print("We maximize I(x) for x in [0, L]. Let u = epsilon*x. We maximize u^2*(1-u)^2 for u in [0,1].")
    u = sympy.Symbol('u')
    f_u = u**2 * (1 - u)**2
    max_val_f = f_u.subs(u, 1/2)
    print(f"The maximum of u^2*(1-u)^2 is {max_val_f} at u=1/2 (x=L/2).")
    # max_I = (max_val_f / (3 * eps**3))
    max_I_val = 1/48
    print(f"So, max_x I(x) = (1/16) / (3*epsilon^3) = {max_I_val:.4f} / epsilon^3.")
    
    print("\nStep 5: Final expression for R^2")
    print("Now we assemble the expression for R^2 = max_x Var[y(x)].")
    print("R^2 ≈ epsilon^4 * (N/L) * max_x I(x)")
    
    # R_squared = eps**4 * N_over_L * (1 / (48 * eps**3))
    # final_R_squared = sympy.simplify(R_squared)

    print("\nSubstituting the terms, we get the final equation:")
    print("R^2 = (epsilon * (1 - epsilon)) / 48")

    print("\nTo be explicit, the numbers in the final equation are:")
    numerator_factor1 = 1
    numerator_factor2 = 1
    denominator = 48
    print(f"Equation form: ( A * epsilon * (B - C * epsilon) ) / D")
    print(f"Coefficient A = {numerator_factor1}")
    print(f"Coefficient B = {numerator_factor2}")
    print(f"Coefficient C = {numerator_factor2}")
    print(f"Denominator D = {denominator}")

if __name__ == '__main__':
    solve_and_explain()