import numpy as np

def analyze_statement_E():
    """
    Analyzes the statement E: "Any strictly convex function has a unique global minimizer."
    This function provides a counterexample to prove the statement is false.
    """
    print("Analyzing Statement E: 'Any strictly convex function has a unique global minimizer.'")
    print("-" * 75)

    # We will use f(x) = e^x as a counterexample.
    # A function is strictly convex if its second derivative is always positive.
    # For f(x) = e^x:
    # First derivative, f'(x) = e^x
    # Second derivative, f''(x) = e^x
    print("Let's consider the function f(x) = e^x.")
    print("The second derivative of f(x) is also e^x, which is greater than 0 for all real numbers x.")
    print("Therefore, f(x) = e^x is a strictly convex function.")
    print("\nNow, let's test if this function has a global minimizer.")
    print("A global minimizer x* would be a point where f(x*) <= f(x) for all other x.")
    
    # We evaluate the function for a sequence of numbers decreasing towards negative infinity.
    test_values = [-1.0, -10.0, -50.0, -100.0, -500.0]
    
    print("\nEvaluating f(x) for a series of decreasing x values:")
    for x in test_values:
        # np.exp can handle the calculation for very small results
        result = np.exp(x)
        print(f"For x = {x: <6}, f(x) = e^({x}) = {result}")
        
    print("\n" + "-" * 75)
    print("As the output shows, as x approaches negative infinity, the value of f(x) gets")
    print("closer and closer to 0, but it never reaches 0 or a negative value.")
    print("There is no single smallest value for the function; we can always find a smaller")
    print("x that gives a smaller f(x). The function's infimum (greatest lower bound) is 0,")
    print("but this value is never attained.")

    print("\nCONCLUSION: The strictly convex function f(x) = e^x does not have a global minimizer.")
    print("This proves that statement E is FALSE.")


if __name__ == '__main__':
    analyze_statement_E()