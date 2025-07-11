import numpy as np
import itertools
from sklearn.linear_model import LogisticRegression

def solve():
    """
    Analyzes which logical operators are learnable by a logistic regression model
    on a given heuristic representation of two embeddings.
    """
    # Let's assume 2D embeddings for simplicity to test mixing-dimension operators.
    # h1 = [p1, q1], h2 = [p2, q2]
    # There are 2^4 = 16 possible input combinations for (p1, q1, p2, q2).
    
    input_configs = list(itertools.product([0, 1], repeat=4))
    
    # Heuristic: h = [h1, h2, |h1-h2|, h1⊙h2]
    # For 2D vectors, this expands to an 8-dimensional feature vector.
    X_features = []
    for p1, q1, p2, q2 in input_configs:
        h1 = np.array([p1, q1])
        h2 = np.array([p2, q2])
        h_diff = np.abs(h1 - h2)
        h_prod = h1 * h2
        heuristic_vector = np.concatenate([h1, h2, h_diff, h_prod])
        X_features.append(heuristic_vector)
    
    X_features = np.array(X_features)

    # Define the logical operators to be tested.
    # We use lambda functions to generate the target labels (y) for each operator.
    # Note: For element-wise ops, we test on the first dimension without loss of generality.
    operators = {
        'X': lambda p1, q1, p2, q2: p1 ^ p2,
        'C': lambda p1, q1, p2, q2: p1 & p2,
        'D': lambda p1, q1, p2, q2: p1 | p2,
        'E': lambda p1, q1, p2, q2: 1 if p1 == p2 else 0,
        'I': lambda p1, q1, p2, q2: 1 if not p1 or p2 else 0,
        'X\'': lambda p1, q1, p2, q2: p1 ^ q2,
        'C\'': lambda p1, q1, p2, q2: p1 & q2,
        'D\'': lambda p1, q1, p2, q2: p1 | q2,
        'E\'': lambda p1, q1, p2, q2: 1 if p1 == q2 else 0,
        'I\'': lambda p1, q1, p2, q2: 1 if not p1 or q2 else 0,
    }

    unlearnable_ops = []

    print("Analyzing learnability of each operator:")
    print("-" * 40)

    for name, op_func in operators.items():
        # Generate target vector y for the current operator
        y_target = np.array([op_func(p1, q1, p2, q2) for p1, q1, p2, q2 in input_configs])

        # If all labels are the same, it's trivially learnable by a constant model.
        if len(np.unique(y_target)) < 2:
            learnable = True
        else:
            # Use a logistic regression model to check for linear separability.
            # A high C value minimizes regularization to find a separator if one exists.
            model = LogisticRegression(solver='liblinear', C=1e9, random_state=42)
            model.fit(X_features, y_target)
            
            # If the accuracy is 1.0, the data is linearly separable and the operator is learnable.
            score = model.score(X_features, y_target)
            learnable = score == 1.0

        print(f"Operator '{name}': {'Learnable' if learnable else 'Not Learnable'} (Accuracy: {score:.2f})")

        if not learnable:
            unlearnable_ops.append(name.replace("'", "'"))
            
    # Sort for consistent output, matching answer format I
    unlearnable_ops.sort()
    
    print("-" * 40)
    print("The following operators cannot be learned with the given heuristic and a linear model:")
    # We join them into a single string to match the format of the answer choices
    final_result = "".join(unlearnable_ops)
    print(final_result)


solve()
<<<I>>>