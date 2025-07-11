import math

def demonstrate_counterexample():
    """
    This function demonstrates that the statement 'Any strictly convex function
    has a unique global minimizer' is false by providing a counterexample.

    The function f(x) = e^x is strictly convex, but it does not have a global minimum
    on the set of real numbers.
    """
    print("Analyzing the strictly convex function f(x) = e^x.")
    print("We will show that it has no minimum value by checking its value for progressively smaller x.")
    print("-" * 40)

    # Let's check decreasing values for x and see the output of f(x)
    test_values = [-1, -10, -100, -1000]

    for x in test_values:
        # Calculate f(x) = e^x
        value = math.exp(x)
        print(f"For x = {x: <5}, f(x) = e^({x}) = {value}")

    print("-" * 40)
    print("As x approaches negative infinity, the value of f(x) approaches 0 but never reaches it.")
    print("Therefore, there is no single value of x for which the function is minimized.")
    print("\nThis demonstrates that the statement 'Any strictly convex function has a unique global minimizer' is FALSE.")


if __name__ == "__main__":
    demonstrate_counterexample()