import numpy as np

def explain_false_statement():
    """
    This script explains why statement E is false by providing a counterexample.
    """
    print("Analyzing the statement: 'E. Any strictly convex function has a unique global minimizer'")
    print("-" * 70)

    print("This statement is FALSE.\n")

    print("Explanation:")
    print("A function is 'strictly convex' if the line segment connecting any two points on its graph lies strictly above the graph.")
    print("While it is true that IF a strictly convex function has a global minimum, that minimum MUST be unique.")
    print("However, a strictly convex function is NOT guaranteed to have a global minimum at all.\n")

    print("Let's consider a counterexample: the function f(x) = e^x (the exponential function).\n")

    # 1. Show strict convexity
    print("1. Is f(x) = e^x strictly convex?")
    print("   Yes. A function is strictly convex if its second derivative is always positive.")
    print("   The first derivative of e^x is e^x.")
    print("   The second derivative of e^x is also e^x.")
    print("   Since e^x > 0 for all real numbers x, the function f(x) = e^x is strictly convex everywhere.\n")

    # 2. Show it has no global minimizer
    print("2. Does f(x) = e^x have a global minimizer?")
    print("   No. A global minimizer is a point 'x' where f(x) is at its lowest possible value.")
    print("   The range of f(x) = e^x is (0, +inf).")
    print("   This means the function's value is always greater than 0.")
    print("   As x approaches negative infinity, f(x) gets closer and closer to 0, but it never actually reaches 0.")
    print("   For any value e^x you pick, I can always find a smaller value (e.g., by choosing a smaller x).")
    print("   Therefore, there is no single point 'x' that produces a minimum value. The function has an infimum (greatest lower bound) of 0, but no minimum.\n")

    print("Conclusion:")
    print("Since f(x) = e^x is a strictly convex function that does not have a global minimizer, the statement 'Any strictly convex function has a unique global minimizer' is proven to be false.")
    print("-" * 70)

if __name__ == '__main__':
    explain_false_statement()