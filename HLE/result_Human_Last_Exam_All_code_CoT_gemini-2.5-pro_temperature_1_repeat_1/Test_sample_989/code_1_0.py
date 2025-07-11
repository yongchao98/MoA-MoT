import numpy as np

def analyze_svm_statements():
    """
    Analyzes several statements about Support Vector Machines (SVMs)
    and identifies the one that is not true.
    """
    print("Analyzing Statement E: 'Any strictly convex function has a unique global minimizer.'")
    print("-" * 70)
    print("This statement is FALSE.")
    print("\nA function is strictly convex if its graph curves upwards. A key property is that IF a minimizer exists, it must be unique.")
    print("However, the existence of a minimizer is not guaranteed.")
    print("\nLet's use the function f(x) = e^x as a counterexample.")
    print("The second derivative of this function is f''(x) = e^x, which is always positive. Therefore, f(x) is strictly convex.")
    print("Now, let's see if it has a global minimum by evaluating it at progressively smaller values of x:")

    x_values = [5, 0, -5, -10, -20, -50, -100]
    print("\nx          f(x) = e^x")
    print("-------------------------")
    for x in x_values:
        # np.exp(x) calculates e to the power of x
        result = np.exp(x)
        print(f"{x:<10} {result:.10e}")

    print("\nAs x approaches negative infinity, f(x) gets closer and closer to 0 but never reaches it.")
    print("There is no single 'x' value that gives the smallest possible 'f(x)'.")
    print("Thus, f(x) = e^x is a strictly convex function that does not have a global minimizer, which proves the statement is false.")
    print("-" * 70)

    print("\nFor completeness, here is why the other statements are TRUE:")
    print("\nA. Mathematically, you cannot have a valid SVM solution using support vectors from only one class.")
    print("   TRUE: This would violate a KKT condition (Σ α_i * y_i = 0) of the SVM dual problem.")

    print("\nB. With imbalanced or asymmetric data, having unequal margins can be optimal for SVM.")
    print("   TRUE: By assigning different misclassification penalties (C) to each class, the SVM can create a decision boundary closer to the majority class to better protect the minority class.")

    print("\nC. Effective mapping to an infinite-dimensional space is computationally tractable for some kernels.")
    print("   TRUE: This is the purpose of the 'kernel trick'. The RBF kernel, for instance, corresponds to an infinite-dimensional feature space, yet its value is easily computed.")

    print("\nD. It is possible to add or move data points and not affect the decision boundary at all, as long as they're interior points.")
    print("   TRUE: The decision boundary is determined solely by the support vectors. Other points (interior points) that are correctly classified and lie outside the margin do not influence it.")


analyze_svm_statements()
<<<E>>>