import sympy

def find_transformation_on_x():
    """
    This script outlines the derivation and provides the general form for the
    Lie group transformation on the spatial variable 'x' for the PDE:
    u_t = u_xx + (k1*ln(u) + k2)*u
    """
    # Define symbolic variables for the problem
    t, x, epsilon = sympy.symbols('t x epsilon', real=True)
    k1 = sympy.Symbol('k1', real=True, nonzero=True)
    c2, c3, c4 = sympy.symbols('c2 c3 c4', real=True)

    print("Step 1: Infinitesimal Generators")
    print("---------------------------------")
    print("Solving the determining equations for the PDE gives the following general forms for the infinitesimals (assuming k1 != 0):")

    # General form of the infinitesimals derived from the symmetry analysis
    tau = c2
    xi_t = c4 - (2 * c3 / k1) * sympy.exp(k1 * t)

    print(f"τ(t) = {tau}")
    print(f"ξ(t) = {xi_t}")
    print("\nHere, c2, c3, and c4 are arbitrary constants.")

    print("\nStep 2: Calculating the Finite Transformation for x")
    print("--------------------------------------------------")
    print("To find the finite transformation x'(ε), we solve the system of ODEs:")
    print("  dt'/dε = τ(t')  with t'(0) = t")
    print("  dx'/dε = ξ(t')  with x'(0) = x")

    # Solve for t'(ε)
    # The equation is dt'/dε = c2, which gives t' = t + c2*ε
    t_prime_expr = t + c2 * epsilon

    print(f"\nSolving for t' gives: t'(ε) = {t_prime_expr}")

    # Substitute t' into ξ to get ξ(t')
    xi_t_prime = xi_t.subs(t, t_prime_expr)

    print(f"The ODE for x' becomes: dx'/dε = {xi_t_prime}")

    # Integrate ξ(t'(ε)) with respect to ε from 0 to ε to find the change in x
    # This integration has two cases based on whether c2 is zero.
    
    # We compute the general integral which is valid for c2 != 0.
    s = sympy.Symbol('s')
    integrand = xi_t.subs(t, t + c2*s)
    integral_part = sympy.integrate(integrand, (s, 0, epsilon))
    
    # The full solution for x'
    x_prime_sol = x + integral_part
    
    print("\nStep 3: General Representation of the Transformation on x")
    print("--------------------------------------------------------")
    print("The final form of the transformation x'(ε) for the general case where c2 is not necessarily zero is:")
    
    # The print function in sympy is used for a clearer math output
    # We build the equation string manually to meet the user's request.
    final_equation = f"x' = x + ({c4})*epsilon - (2*{c3}/(({k1}**2)*{c2})) * exp({k1}*t) * (exp({k1}*{c2}*epsilon) - 1)"
    
    print("\nFinal Equation:")
    sympy.pprint(sympy.Eq(sympy.Symbol("x'"), x_prime_sol), use_unicode=False)
    
    print("\nThe above equation with all parameters explicitly shown is:")
    print(final_equation)
    print("\nNote: This general form correctly simplifies to the case for c2=0 by taking the limit, where it becomes x' = x + ε * ξ(t).")

if __name__ == '__main__':
    find_transformation_on_x()