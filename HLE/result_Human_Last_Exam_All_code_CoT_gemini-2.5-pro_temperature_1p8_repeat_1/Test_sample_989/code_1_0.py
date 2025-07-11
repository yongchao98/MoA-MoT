import numpy as np

def explain_false_statement():
    """
    This script explains why one of the statements about Support Vector Machines is false.
    The goal is to analyze the following statement:

    E. Any strictly convex function has a unique global minimizer.

    This statement is a general mathematical claim, and we can prove it false with a counterexample.
    """

    print("Analyzing Statement E: 'Any strictly convex function has a unique global minimizer.'")
    print("-" * 70)

    # Explanation
    print("This statement is FALSE.\n")
    print("Reasoning: While a strictly convex function can have at most one global minimizer (if one exists, it is unique), its existence is not guaranteed.")
    print("\nCounterexample: The function f(x) = e^x (exponential function).")

    # 1. Prove strict convexity
    print("\nStep 1: Show that f(x) = e^x is strictly convex.")
    print("A function is strictly convex if its second derivative is always positive.")
    print("f'(x)  = d/dx(e^x) = e^x")
    print("f''(x) = d/dx(e^x) = e^x")
    print("Since e^x is always greater than 0 for any real number x, the function f(x) = e^x is strictly convex everywhere.")

    # 2. Show it has no minimum
    print("\nStep 2: Show that f(x) = e^x has no global minimum.")
    print("The function e^x is always positive.")
    print("As x approaches negative infinity, f(x) approaches 0, but it never actually reaches 0.")
    print(f"For example, f(-10) = {np.exp(-10):.5f}")
    print(f"For example, f(-100) = {np.exp(-100):.5f}")
    print("We can always find a smaller value for x that results in a smaller value for f(x).")
    print("Therefore, there is no single value of x for which f(x) is the minimum. The function has an infimum (greatest lower bound) of 0, but it never attains this value.")

    # Conclusion
    print("\n" + "-" * 70)
    print("Conclusion: Since we have found a counterexample (f(x) = e^x), the general statement that 'Any strictly convex function has a unique global minimizer' is not true.")

if __name__ == '__main__':
    explain_false_statement()