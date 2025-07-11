import numpy as np

def demonstrate_false_statement():
    """
    This function demonstrates a counterexample to the statement:
    'Any strictly convex function has a unique global minimizer'.

    The function f(x) = e^x is strictly convex, but it does not have a global minimum
    on the set of real numbers. Its greatest lower bound (infimum) is 0, but it never
    actually reaches this value for any finite x.

    We will show this by evaluating f(x) for a series of large negative numbers.
    """
    print("Analyzing Statement E: 'Any strictly convex function has a unique global minimizer'")
    print("This statement is FALSE. Here is a demonstration with a counterexample.\n")
    print("Consider the function f(x) = e^x, which is strictly convex.")
    print("If it had a global minimum, there would be a value x* where f(x*) is the smallest possible value.")
    print("Let's evaluate the function as x gets more negative:\n")

    x_values = [-1, -10, -50, -100]
    for x in x_values:
        # The final equation is y = e^x
        y = np.exp(x)
        # Output each number in the final equation as requested
        print(f"For x = {x}, the function value is f({x}) = e^({x}) = {y}")

    print("\nAs you can see, as x decreases, the function value f(x) gets closer to 0.")
    print("However, there is no finite value of x for which f(x) actually equals 0.")
    print("Therefore, the strictly convex function f(x) = e^x does not have a global minimum.")
    print("This proves the statement is false.")

demonstrate_false_statement()