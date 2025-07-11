import numpy as np

def demonstrate_false_statement():
    """
    This function demonstrates that statement E is false.
    Statement E: Any strictly convex function has a unique global minimizer.
    We use the function f(x) = e^x as a counterexample.
    """

    # Define the function and its second derivative
    f = lambda x: np.exp(x)
    f_second_derivative = lambda x: np.exp(x)

    print("--- Demonstrating Statement E is False ---")
    print("Statement E: 'Any strictly convex function has a unique global minimizer.'\n")
    print("Let's consider the function f(x) = e^x.")
    print("A function is strictly convex if its second derivative is always strictly positive (> 0).")
    print("The second derivative of f(x) is also e^x, which is > 0 for all real numbers x.")
    print("Therefore, f(x) = e^x is a strictly convex function.\n")
    print("Now, let's see if it has a global minimum by observing its value as x becomes more negative:")

    # Show that the function value continuously decreases as x -> -infinity
    x_values = [10, 0, -10, -20, -50, -100]
    for x in x_values:
        print(f"For x = {x:5.0f}, f(x) = e^x = {f(x):.10f}")

    print("\nAs x approaches negative infinity, f(x) approaches 0 but never reaches it.")
    print("This means the function has no single point at which it achieves a minimum value.")
    print("It has an infimum (greatest lower bound) of 0, but no minimum.")
    print("\nConclusion: We have found a strictly convex function (e^x) that does not have a global minimizer.")
    print("Therefore, the statement 'Any strictly convex function has a unique global minimizer' is not true.")


if __name__ == '__main__':
    demonstrate_false_statement()
