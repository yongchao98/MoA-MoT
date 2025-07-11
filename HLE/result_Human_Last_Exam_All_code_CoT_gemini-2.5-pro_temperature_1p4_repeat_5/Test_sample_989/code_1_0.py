import numpy as np

def analyze_svm_statements():
    """
    Analyzes each statement about Support Vector Machines to find the false one.
    """
    print("Analyzing the statements to find the one that is not true:\n")

    print("--- Statement A ---")
    print("Statement: Mathematically, you cannot have a valid SVM solution using support vectors from only one class.")
    print("Analysis: This is TRUE. The SVM dual formulation has a constraint `sum(alpha_i * y_i) = 0`, where `alpha_i > 0` only for support vectors and `y_i` is the class label (+1 or -1). For this sum to be zero, there must be positive `alpha_i` associated with both positive and negative `y_i`. Therefore, support vectors must come from both classes.")
    print("-" * 20 + "\n")

    print("--- Statement B ---")
    print("Statement: With imbalanced or asymmetric data, having unequal margins can be optimal for SVM.")
    print("Analysis: This is TRUE. In class-weighted SVMs, different misclassification penalties (C_positive, C_negative) are used for each class. This leads the model to create an asymmetric margin that is wider on the majority class side and narrower on the minority class side to reduce misclassifications of the minority class.")
    print("-" * 20 + "\n")

    print("--- Statement C ---")
    print("Statement: Effective mapping to an infinite-dimensional space is computationally tractable for some kernels.")
    print("Analysis: This is TRUE. This is the main benefit of the 'kernel trick'. The Radial Basis Function (RBF) kernel, K(x, z) = exp(-gamma * ||x-z||^2), corresponds to a mapping into an infinite-dimensional feature space, but we never compute the mapping itself, only the kernel function, which is computationally cheap.")
    print("-" * 20 + "\n")

    print("--- Statement D ---")
    print("Statement: It is possible to add or move data points and not affect the decision boundary at all, as long as they're interior points.")
    print("Analysis: This is TRUE. The decision boundary is determined exclusively by the support vectors (the points on or inside the margin). Points correctly classified and far from the boundary (interior points) have no impact on the hyperplane's location.")
    print("-" * 20 + "\n")

    print("--- Statement E ---")
    print("Statement: Any strictly convex function has a unique global minimizer.")
    print("Analysis: This is FALSE. While a strictly convex function can have at most one global minimizer, it is not guaranteed to have one at all. A classic counterexample is the function f(x) = e^x.")
    print("  - Convexity Check: The second derivative is f''(x) = e^x, which is always positive. This means the function is strictly convex.")
    print("  - Minimizer Check: The function approaches 0 as x approaches negative infinity (infimum is 0), but it never reaches a minimum value. For any point x, we can always find a point x-1 where the function value is smaller (e^(x-1) < e^x).")
    print("-" * 20 + "\n")

    print("Conclusion: Statement E is the only one that is not true.")

if __name__ == '__main__':
    analyze_svm_statements()
    print("<<<E>>>")