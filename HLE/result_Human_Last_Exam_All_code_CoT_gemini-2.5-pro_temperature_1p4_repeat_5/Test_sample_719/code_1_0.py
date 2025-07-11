import sympy as sp

def solve_linearized_flow():
    """
    This function uses symbolic mathematics to derive the expression for theta'(t).
    It follows the steps outlined in the explanation.
    """
    # 1. Define all necessary symbolic functions and variables.
    t = sp.Symbol('t')
    f = sp.Function('f')(t)
    c1 = sp.Function('c1')(t)
    c2 = sp.Function('c2')(t)
    theta = sp.Function('theta')(t)
    r = sp.Function('r')(t)

    print("Step 1: Define the system of ODEs for coefficients c1(t) and c2(t).")
    # 2. Define the system of ODEs for c1 and c2 based on the derivation.
    # c1' = -c1 * f' / f
    # c2' = c1 * f
    f_prime = f.diff(t)
    c1_prime_expr = -c1 * f_prime / f
    c2_prime_expr = c1 * f
    print(f"c1'(t) = {c1_prime_expr}")
    print(f"c2'(t) = {c2_prime_expr}\n")

    print("Step 2: Define theta'(t) in terms of c1(t), c2(t) and their derivatives.")
    # 3. Use the formula for the derivative of arctan(c2/c1).
    # theta' = (c1*c2' - c2*c1') / (c1^2 + c2^2)
    theta_prime_formula = (c1 * c2.diff(t) - c2 * c1.diff(t)) / (c1**2 + c2**2)
    print(f"theta'(t) = {theta_prime_formula}\n")
    
    print("Step 3: Substitute the ODEs into the formula for theta'(t).")
    # 4. Substitute the expressions for the derivatives.
    theta_prime_intermediate = theta_prime_formula.subs({
        c1.diff(t): c1_prime_expr,
        c2.diff(t): c2_prime_expr
    })
    print(f"theta'(t) = {sp.simplify(theta_prime_intermediate)}\n")
    
    print("Step 4: Substitute polar coordinate expressions c1=r*cos(theta), c2=r*sin(theta).")
    # 5. Substitute the polar coordinate mapping and simplify.
    final_expr = theta_prime_intermediate.subs({
        c1: r * sp.cos(theta),
        c2: r * sp.sin(theta)
    })
    final_expr_simplified = sp.simplify(final_expr)

    print("--- Final Result ---")
    print("The expression for theta'(t) is:")
    sp.pprint(final_expr_simplified)
    print("\nThis can be written as:")
    print("theta'(t) = f(t)*cos^2(theta(t)) + (f'(t)/f(t))*cos(theta(t))*sin(theta(t))")

    # The user asked to output each number in the final equation.
    # Since the equation is symbolic, we will describe its components.
    print("\nBreaking down the final equation:")
    # Using Poly to extract coefficients, a bit overkill but robust
    cos2_term = f * sp.cos(theta)**2
    sincos_term = f_prime/f * sp.cos(theta) * sp.sin(theta)
    
    print(f"The term with cos^2(theta(t)) has a coefficient of: f(t)")
    print(f"The term with cos(theta(t))*sin(theta(t)) has a coefficient of: f'(t)/f(t)")
    print("All other terms are zero.")


if __name__ == '__main__':
    solve_linearized_flow()