import sympy

def solve_fixed_point():
    """
    This function derives the leading order expression for the fixed point coupling u*
    in phi^4 theory near d=4 dimensions.
    """
    # Define the necessary symbols for our calculation
    epsilon, pi = sympy.symbols('epsilon pi', real=True, positive=True)
    lamda = sympy.symbols('lambda')

    # Step 1: Define the one-loop beta function for the coupling λ.
    # The interaction term in the Lagrangian is (λ/4!)ϕ⁴.
    # The beta function in d = 4 - ϵ dimensions is β(λ) = -ϵλ + (3/(16π²))λ².
    beta_lam_expr = -epsilon * lamda + (3 / (16 * pi**2)) * lamda**2

    # Step 2: Solve for the fixed point λ* by setting β(λ*) = 0.
    # We are interested in the non-trivial (Wilson-Fisher) fixed point.
    fixed_point_solutions = sympy.solve(beta_lam_expr, lamda)
    lamda_star = fixed_point_solutions[1]  # The non-zero solution is at index 1

    # Step 3: Relate λ to the coupling u.
    # The coupling 'u' in an interaction term uϕ⁴ corresponds to λ/4!
    # 4! = 24
    u_star = lamda_star / 24
    
    # Step 4: Print the derivation in a clear, step-by-step format.
    print("Derivation of the fixed point coupling u*:")
    print("-" * 50)
    
    print("The beta function for the coupling λ from the (λ/4!)ϕ⁴ term is:")
    print(f"  β(λ) = -ϵ*λ + (3/(16*π²))*λ²")
    
    print("\nSolving β(λ*) = 0 gives the Wilson-Fisher fixed point λ*:")
    # We represent the equation by showing the non-trivial result
    print(f"  λ* = {sympy.pretty(lamda_star, use_unicode=False)}")

    print("\nThe coupling u (from uϕ⁴) is related to λ by u = λ / 4! = λ / 24.")
    print("Substituting λ* gives the fixed point u*:")
    print(f"  u* = λ* / 24 = ({sympy.pretty(lamda_star, use_unicode=False)}) / 24")
    
    # Simplify the final expression for u*
    u_star_simplified = sympy.simplify(u_star)
    
    # As requested, output the numbers in the final equation u* = C * π² * ϵ
    # C = 2/9
    num, den = sympy.fraction(u_star_simplified / (pi**2 * epsilon))

    print("\nFinally, simplifying the expression gives the result:")
    final_equation = f"  u* = ({num} * π² / {den}) * ϵ"
    print(final_equation)
    print("\nEach number in the final equation is:")
    print(f"  Numerator of the coefficient: {num}")
    print(f"  Denominator of the coefficient: {den}")
    print("  The expression involves π to the power of: 2")
    print("  The expression is linear in ϵ (power: 1)")

if __name__ == '__main__':
    solve_fixed_point()