import numpy as np

def analyze_svm_statements():
    """
    Provides a detailed analysis of each statement about Support Vector Machines (SVMs)
    and demonstrates why one of them is false.
    """
    print("Analyzing the provided statements about SVMs:")
    print("="*45)

    print("\nA. Mathematically, you cannot have a valid SVM solution using support vectors from only one class.")
    print("   Analysis: TRUE. The decision boundary is placed to maximize the margin between classes, which requires 'pushback' from points in both classes. These points are the support vectors.")

    print("\nB. With imbalanced or asymmetric data, having unequal margins can be optimal for SVM.")
    print("   Analysis: TRUE. By using different penalty parameters (C) for each class (cost-sensitive SVM), one can effectively create asymmetric margins to better handle imbalanced data.")

    print("\nC. Effective mapping to an infinite-dimensional space is computationally tractable for some kernels.")
    print("   Analysis: TRUE. This is the 'kernel trick'. The RBF kernel, for example, maps to an infinite-dimensional space, but we only ever need to compute the scalar kernel value, which is efficient.")

    print("\nD. It is possible to add or move data points and not affect the decision boundary at all, as long as they're interior points.")
    print("   Analysis: TRUE. The SVM boundary is determined only by support vectors. Interior points (correctly classified and outside the margin) do not influence the boundary's position.")

    print("\nE. Any strictly convex function has a unique global minimizer.")
    print("   Analysis: FALSE. A strictly convex function is NOT guaranteed to have a global minimum. If a minimum exists, it is unique. The statement is therefore false.")
    
    print("\n--- Demonstration for Statement E ---")
    # A strictly convex function is f(x) = e^x. Its second derivative is e^x, which is > 0 for all x.
    # However, it has no global minimum.
    f = np.exp
    print("Let's consider the strictly convex function f(x) = e^x.")
    print("The function's value approaches 0 as x goes to negative infinity but never reaches it.")
    print(f"f(-10)   = {f(-10):.6f}")
    print(f"f(-100)  = {f(-100):.6e}")
    print(f"f(-1000) = {f(-1000):.6e}")
    print("Since the function never reaches its infimum (0), it has no global minimum.")

    print("\n" + "="*45)
    print("Conclusion: Statement E is the only one that is not true.")

if __name__ == '__main__':
    analyze_svm_statements()