import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import NotFittedError

def analyze_and_identify_unlearnable_operators():
    """
    Analyzes a list of logical operators to determine if they can be learned
    by a logistic regression model given a specific feature heuristic.
    """

    def check_learnability(features, labels, op_name):
        """
        Checks if a logistic regression model can perfectly learn the mapping.
        A perfect score (1.0) means the function is linearly separable in the given feature space.
        """
        # If all labels are the same, the classification is trivial.
        if len(np.unique(labels)) < 2:
            return True
            
        model = LogisticRegression(solver='liblinear', C=1e5, random_state=0)
        model.fit(features, labels)
        score = model.score(features, labels)
        return score == 1.0

    unlearnable_operators = []
    
    # Define the 4 possible input pairs for two boolean variables (p, q)
    # p could be h1[i], q could be h2[i] or h2[j]
    inputs = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    p = inputs[:, 0]
    q = inputs[:, 1]

    # --- Part 1: Analyze Element-wise Operators ---
    print("--- Analyzing Element-wise Operators f(h1[i], h2[i]) ---")
    print("Feature vector: [h1[i], h2[i], |h1[i]-h2[i]|, h1[i]*h2[i]]\n")

    # The 4x4 feature matrix for the 4 possible input pairs
    features_elementwise = np.array([
        p,                      # h1[i]
        q,                      # h2[i]
        np.abs(p - q),          # |h1[i]-h2[i]|
        p * q                   # h1[i]*h2[i]
    ]).T

    elementwise_ops = {
        'X (XOR)': lambda p, q: np.logical_xor(p, q),
        'C (Conjunction)': lambda p, q: np.logical_and(p, q),
        'D (Disjunction)': lambda p, q: np.logical_or(p, q),
        'E (Equivalence)': lambda p, q: np.equal(p, q),
        'I (Implication)': lambda p, q: np.logical_or(np.logical_not(p), q)
    }

    for name, op_func in elementwise_ops.items():
        labels = op_func(p, q).astype(int)
        learnable = check_learnability(features_elementwise, labels, name)
        if not learnable:
            unlearnable_operators.append(name.split(' ')[0])
        print(f"Operator: {name.ljust(18)} | Learnable: {learnable}")

    # --- Part 2: Analyze Mixing-Dimension Operators ---
    print("\n--- Analyzing Mixing-Dimension Operators f(h1[i], h2[j]) ---")
    print("Feature vector: [h1[i], h2[j]] (plus other irrelevant features)\n")

    # The 4x2 feature matrix for the 4 possible input pairs
    features_mixing = inputs

    mixing_ops = {
        "X' (XOR)": lambda p, q: np.logical_xor(p, q),
        "C' (Conjunction)": lambda p, q: np.logical_and(p, q),
        "D' (Disjunction)": lambda p, q: np.logical_or(p, q),
        "E' (Equivalence)": lambda p, q: np.equal(p, q),
        "I' (Implication)": lambda p, q: np.logical_or(np.logical_not(p), q)
    }

    for name, op_func in mixing_ops.items():
        labels = op_func(p, q).astype(int)
        learnable = check_learnability(features_mixing, labels, name)
        if not learnable:
            # Extract the short name like "X'" from "X' (XOR)"
            unlearnable_operators.append(name.split(' ')[0])
        print(f"Operator: {name.ljust(18)} | Learnable: {learnable}")

    # --- Part 3: Final Result ---
    print("\n-----------------------------------------------------------")
    print("Conclusion: The operators that cannot be learned are those that are not linearly")
    print("separable in their respective feature spaces.")
    print("\nFinal list of unlearnable operators:", ''.join(sorted(unlearnable_operators)))
    print("-----------------------------------------------------------")

if __name__ == '__main__':
    analyze_and_identify_unlearnable_operators()