import numpy as np

def explain_false_statement():
    """
    Explains why statement E is false using a mathematical counterexample.
    """
    print("The statement that is not true is:")
    print("E. Any strictly convex function has a unique global minimizer\n")
    print("-" * 70)
    print("EXPLANATION")
    print("-" * 70)
    print("This statement is FALSE because it claims a global minimizer always exists, which is not guaranteed.")
    print("A statement is proven false by finding a single counterexample.\n")

    print("Counterexample: The function f(x) = e^x (numpy.exp(x))\n")

    # Step 1: Prove strict convexity
    print("1. Is the function f(x) = e^x strictly convex?")
    print("   A common way to check for strict convexity is to examine the second derivative.")
    print("   If the second derivative is strictly positive (f''(x) > 0) for all x, the function is strictly convex.")
    print("   - First derivative (f'(x)) of e^x is e^x.")
    print("   - Second derivative (f''(x)) of e^x is also e^x.")
    print("   The function e^x is greater than 0 for any real number x. Therefore, f(x) = e^x is strictly convex.\n")

    # Step 2: Show that a global minimizer does not exist
    print("2. Does the function f(x) = e^x have a global minimizer?")
    print("   A global minimizer is a point x* for which f(x*) is less than or equal to f(x) for all other x.")
    print("   The range of f(x) = e^x is (0, +infinity).")
    print("   As x approaches negative infinity, f(x) gets closer and closer to 0.")
    x_values = np.array([-1, -10, -100, -1000])
    y_values = np.exp(x_values)
    for x, y in zip(x_values, y_values):
        print(f"   f({x}) = {y:.2e}")
    print("\n   The function approaches a lower bound of 0, but it never reaches it.")
    print("   There is no real number x* such that e^(x*) = 0. The infimum is 0, but it is not a minimum.\n")

    # Step 3: Conclusion
    print("Conclusion:")
    print("f(x) = e^x is a strictly convex function that does not have a global minimizer.")
    print("This proves the general statement 'Any strictly convex function has a unique global minimizer' is FALSE.")

if __name__ == "__main__":
    explain_false_statement()
