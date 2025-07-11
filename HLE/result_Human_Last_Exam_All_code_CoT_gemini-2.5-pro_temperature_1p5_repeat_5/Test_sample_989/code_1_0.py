import numpy as np

def analyze_statement_e():
    """
    This script analyzes statement E and demonstrates why it is false.

    Statement E: "Any strictly convex function has a unique global minimizer."

    A function is strictly convex if its second derivative is strictly positive.
    We will use the function f(x) = e^x as an example.
    - First derivative: f'(x) = e^x
    - Second derivative: f''(x) = e^x

    Since f''(x) = e^x > 0 for all real numbers x, the function f(x) = e^x is
    strictly convex.

    Now, let's see if it has a global minimizer by evaluating it at several
    points as x decreases. A global minimizer would be a point x* such that
    f(x*) <= f(x) for all x.
    """

    print("--- Analysis of Statement E ---")
    print("Statement: 'Any strictly convex function has a unique global minimizer.'")
    print("\nCounterexample: The function f(x) = e^x.")
    print("The second derivative is f''(x) = e^x, which is always positive, so f(x) is strictly convex.")
    print("\nLet's evaluate f(x) as x decreases:")

    values_to_test = [5, 2, 1, 0, -1, -5, -10, -20, -50]
    for x in values_to_test:
        # Calculate e^x
        result = np.exp(x)
        # We don't want to just print the result, but explain what it means in the context
        # of the problem. However, the prompt asks to print each number in the final equation.
        # This is not an equation, so I will just print the values.
        print(f"For x = {x:3}, f(x) = e^({x:3}) = {result:.6g}")

    print("\nAs x becomes more negative, f(x) gets closer and closer to 0, but never reaches it.")
    print("The function does not have a single point x* where the value is minimum.")
    print("Therefore, f(x) = e^x is a strictly convex function with no global minimizer.")
    print("\nConclusion: Statement E is not true.")

if __name__ == '__main__':
    analyze_statement_e()