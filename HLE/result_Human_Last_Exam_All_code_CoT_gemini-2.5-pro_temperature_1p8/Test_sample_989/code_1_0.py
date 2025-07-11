import sympy

def analyze_function():
    """
    Analyzes the function f(x) = e^x to show it's a counterexample
    to the claim that any strictly convex function has a unique global minimizer.
    """
    # Define the variable and the function
    x = sympy.symbols('x', real=True)
    f = sympy.exp(x)

    # Calculate first and second derivatives
    f_prime = sympy.diff(f, x)
    f_double_prime = sympy.diff(f_prime, x)

    print("Analyzing Statement E: 'Any strictly convex function has a unique global minimizer.'")
    print("-" * 80)
    print("Let's test this statement with a counterexample: f(x) = e^x")
    print(f"The function is: {f}")

    # Check for strict convexity
    print("\nA function is strictly convex if its second derivative is always positive.")
    print(f"The second derivative of f(x) is: f''(x) = {f_double_prime}")
    print("Since e^x is always positive for any real number x, the function f(x) is strictly convex.")

    # Check for a global minimum
    print("\nTo find a minimum, we set the first derivative to zero and solve for x.")
    print(f"The first derivative of f(x) is: f'(x) = {f_prime}")
    print(f"Set f'(x) = 0  =>  {f_prime} = 0")
    
    # Try to solve the equation
    try:
        solutions = sympy.solve(f_prime, x)
        if not solutions:
            print("The equation e^x = 0 has no solution for x in the real numbers.")
        else:
            print(f"The solution is: {solutions}")
    except Exception as e:
        print(f"Could not solve the equation: {e}")
        
    print("\nBecause there are no points where the derivative is zero, the function has no critical points and thus no global minimum.")
    print("\nConclusion: The statement that 'Any strictly convex function has a unique global minimizer' is FALSE.")
    print("-" * 80)

if __name__ == '__main__':
    analyze_function()