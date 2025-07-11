import numpy as np

def demonstrate_false_statement():
    """
    This function demonstrates that statement E is false by providing a counterexample.
    
    Statement E: Any strictly convex function has a unique global minimizer.
    
    We will analyze the function f(x) = e^x.
    A function is strictly convex if its second derivative is always positive.
    The first derivative of f(x) is f'(x) = e^x.
    The second derivative of f(x) is f''(x) = e^x.
    
    Since e^x is always greater than 0 for any real number x, f(x) = e^x is a strictly convex function.
    
    Now, let's see if it has a global minimizer. A global minimizer is a point x* where f(x*) is less than or equal to f(x) for all x.
    We will evaluate the function for a sequence of decreasing x values.
    """

    print("--- Analyzing Statement E: 'Any strictly convex function has a unique global minimizer.' ---")
    print("Let's use the function f(x) = e^x as a counterexample.")
    print("The second derivative of f(x) is e^x, which is always > 0. Thus, f(x) is strictly convex.")
    print("\nNow let's evaluate f(x) for some decreasing values of x:")

    test_values = [-1, -10, -50, -100, -1000]
    
    for x in test_values:
        # Calculate f(x) = e^x
        y = np.exp(x)
        print(f"For x = {x:6}, f(x) = e^({x:6}) = {y:.4e}")

    print("\nAs you can see, as x decreases towards negative infinity, f(x) gets closer and closer to 0.")
    print("However, there is no real number x for which f(x) is actually 0 or a negative number.")
    print("The function approaches an infimum (greatest lower bound) of 0, but it never reaches a minimum value.")
    print("\nConclusion: We have found a strictly convex function (e^x) that does not have a global minimizer.")
    print("Therefore, the statement 'Any strictly convex function has a unique global minimizer' is FALSE.")

if __name__ == '__main__':
    demonstrate_false_statement()