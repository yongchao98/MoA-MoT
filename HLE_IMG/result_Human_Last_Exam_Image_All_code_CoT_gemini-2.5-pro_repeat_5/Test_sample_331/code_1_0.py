import sympy as sp
import numpy as np

def solve_problem():
    """
    Solves the multi-step physics and mathematics problem.
    """
    # Step 1: Correspondence between geometries and Hamiltonians
    # H1 ~ (q^3-3p^2q)^2 has D6 symmetry -> F (Hexagon)
    # H2 ~ -p^2q^2 has 4 saddles at (±1, ±1) -> D (Diamond)
    # H3 ~ q^3-3p^2q has D3 symmetry -> C (Triangle)
    # H4 ~ p^2-q^4/4 has 2 saddles on q-axis -> B (Lemon)
    # H5 ~ -cos(4theta) potential has D4 symmetry -> E (Square)
    # H6 ~ -q^3 term breaks up-down symmetry -> A (Teardrop)
    
    n = {'A': 6, 'B': 4, 'C': 3, 'D': 2, 'E': 5, 'F': 1}
    print("Step 1: Correspondence")
    for key, val in n.items():
        print(f"  n_{key} = {val}")
    print("-" * 20)

    # Step 2: Calculate parameters
    n_A, n_B, n_C, n_D, n_E, n_F = n['A'], n['B'], n['C'], n['D'], n['E'], n['F']
    
    # x0 for y(x0)=0
    x0_val = sp.Rational(n_F, n_E)
    
    # beta for Caputo derivative
    beta = sp.Rational(n_E, n_B)
    
    # nu for Riemann-Liouville integral
    nu = sp.Rational(n_C, n_A)
    
    # Find n_max: maximizes T_n(1/n_D) / T_n(0)
    # This simplifies to maximizing the product of exponents in the _2F_1 function.
    # H6 (A2 singularity, p^2+q^3) has exponents (1/2, 1/3) -> product 1/6 (Largest)
    n_max = 6
    
    # Find n_S3_min: index of 3rd smallest integral over boundary oscillation disk
    # Order of complexity of singularities: A2 < A3 < D3 < D4 < D5 < D6
    # H6(A2) < H4(A3) < H3(D3) < H2(D4) ...
    # So, 1st is n_A=6, 2nd is n_B=4, 3rd is n_C=3
    # Wait, the indices are H6, H4, H3.
    n_S3_min = 3
    
    print("Step 2: Parameters")
    print(f"  x0 = n_F / n_E = {n_F} / {n_E} = {x0_val}")
    print(f"  lambda = 1 (since r_max_sq is 2 for both H{n_E} and H{n_B})")
    print(f"  nu = n_C / n_A = {n_C} / {n_A} = {nu}")
    print(f"  beta = n_E / n_B = {n_E} / {n_B} = {beta}")
    print(f"  n_max = {n_max}")
    print(f"  n_S3_min = {n_S3_min}")
    print("-" * 20)
    
    # Step 3: Set up and solve the equation for mu
    # The condition y(x0)=0 implies f''(x0)/f'(x0) = g''(x0)/g'(x0)
    
    # 3a: Analyze g(x) = K((lambda*x)^2)^mu = K(x^2)^mu
    # K(alpha) = I^nu T_n_max(alpha). For alpha -> 0, T is constant, so K(alpha) ~ alpha^nu.
    # K(z) ~ z^(1/2). So g(x) ~ (x^2)^(mu/2) = x^mu.
    # This is a key simplifying assumption based on x0 being small.
    x, mu = sp.symbols('x mu')
    # g = x**mu  # Proportionality constant cancels
    # g_prime = sp.diff(g, x)
    # g_prime2 = sp.diff(g_prime, x)
    # g_ratio = sp.simplify(g_prime2 / g_prime)
    # This gives (mu-1)/x
    g_ratio_expr = (mu - 1) / x
    
    # 3b: Analyze f(x)
    # f(x) = D^beta H_{n_S3_min}(n_F, x) = D^(5/4) H3(1, x)
    # H3(p,q) = 1/2 * (-sqrt(1/27)*(q^3 - 3*p**2*q) + p**2 + q**2)
    p_val, q_val = sp.S(n_F), x
    H3 = sp.S(1)/2 * (-sp.sqrt(sp.S(1)/27)*(q_val**3 - 3*p_val**2*q_val) + p_val**2 + q_val**2)
    H3_poly = sp.poly(H3, x)

    # Caputo D^beta f(x) = I^(m-beta) D^m f(x) where m = ceil(beta)
    m = sp.ceiling(beta)
    P_double_prime = sp.diff(H3_poly.as_expr(), x, m)
    
    # f(x) = I^(m-beta) P''(x)
    # Since P'' is a polynomial, we can integrate term-by-term
    # I^alpha(x^k) = Gamma(k+1)/Gamma(k+1+alpha) * x^(k+alpha)
    f_expr = 0
    alpha_I = m - beta
    for term_power, coeff in sp.poly(P_double_prime, x).terms():
        k = term_power[0]
        f_expr += coeff * (sp.gamma(k+1) / sp.gamma(k+1+alpha_I)) * x**(k+alpha_I)

    # Now calculate f'(x) and f''(x)
    f_prime = sp.diff(f_expr, x)
    f_prime2 = sp.diff(f_prime, x)
    
    # The ratio f''/f'
    f_ratio_expr = sp.simplify(f_prime2 / f_prime)

    # 3c: Equate the ratios at x = x0
    g_ratio_at_x0 = g_ratio_expr.subs(x, x0_val)
    f_ratio_at_x0 = f_ratio_expr.subs(x, x0_val)
    
    equation = sp.Eq(g_ratio_at_x0, f_ratio_at_x0)
    
    print("Step 3: The Equation for mu")
    print("  Condition y(x0)=0 leads to f''(x0)/f'(x0) = g''(x0)/g'(x0)")
    print("  Using K(z) ~ z^(1/2), the g-ratio simplifies to (mu-1)/x")
    print(f"  g''/g' at x={x0_val}: {g_ratio_at_x0}")
    print(f"  f''/f' at x={x0_val}: {f_ratio_at_x0.evalf()}")
    print("  Equation to solve:")
    # Pretty print the equation
    sp.pretty_print(equation)
    
    # Solve for mu
    mu_sol = sp.solve(equation, mu)
    
    final_mu = mu_sol[0]
    
    print("\nStep 4: Solution")
    print(f"  The symbolic solution for mu is: {final_mu}")
    print(f"  The numerical value for mu is: {final_mu.evalf()}")
    
    print("\nFinal equation with all numbers:")
    print(f"  ({mu_sol[0].evalf()} - 1) / {x0_val} = {f_ratio_at_x0.evalf()}")
    print(f"  {g_ratio_at_x0.subs(mu, mu_sol[0]).evalf()} = {f_ratio_at_x0.evalf()}")

    return final_mu.evalf()

mu_value = solve_problem()
# The final output needs to be in the specified format
# print(f"\n<<<{mu_value}>>>")
# The problem asks for the answer directly.
# Let's format it to a reasonable number of digits.
final_answer = round(mu_value, 4)
print(f"\nFinal Answer: \n{final_answer}")
print(f'<<<{final_answer}>>>')
