import numpy as np
from sklearn.linear_model import LogisticRegression
import itertools

def get_heuristic_features(h1, h2):
    """Computes the heuristic [h1, h2, |h1-h2|, h1⊙h2]."""
    if not isinstance(h1, np.ndarray):
        h1 = np.array(h1)
    if not isinstance(h2, np.ndarray):
        h2 = np.array(h2)
    
    diff = np.abs(h1 - h2)
    prod = h1 * h2
    return np.concatenate([h1, h2, diff, prod])

def check_learnability(operator_name, target_func, dims):
    """
    Tests if an operator is learnable by checking for linear separability.
    
    Args:
        operator_name (str): The name of the operator (e.g., "C'").
        target_func (callable): A function that takes (h1, h2) and returns the target value.
        dims (int): The number of dimensions for the input embeddings (1 for element-wise, 2 for mixing).

    Returns:
        bool: True if learnable, False otherwise.
    """
    # 1. Generate all possible input embeddings h1, h2
    possible_vectors = list(itertools.product([0, 1], repeat=dims))
    all_h1 = [np.array(v) for v in possible_vectors]
    all_h2 = [np.array(v) for v in possible_vectors]
    
    X_features = []
    y_targets = []

    # 2. Generate dataset of heuristic features and target labels
    for h1 in all_h1:
        for h2 in all_h2:
            # For 1D case, this will be just a single vector pair.
            # For 2D case, this iterates through all 16 h1/h2 vector pairs.
            if dims == 1:
                # Re-wrap h1/h2 to handle scalar vs array logic
                features = get_heuristic_features([h1[0]], [h2[0]])
                target = target_func(h1[0], h2[0])
            else: # dims == 2
                features = get_heuristic_features(h1, h2)
                target = target_func(h1, h2)
            
            X_features.append(features)
            y_targets.append(target)
    
    X = np.array(X_features)
    y = np.array(y_targets)

    # Check if there is more than one class in the output
    if len(np.unique(y)) < 2:
        return True # A function with a constant output is trivially learnable

    # 3. Test for linear separability with a Logistic Regression model
    # We use a very high C to minimize regularization and test for separability
    model = LogisticRegression(solver='liblinear', C=1e10, random_state=0)
    model.fit(X, y)
    
    # If the model gets a perfect score, the data is linearly separable
    score = model.score(X, y)
    return score == 1.0


def main():
    """
    Main function to evaluate all operators and print the unlearnable ones.
    """
    operators = {
        # Element-wise operators (dims=1)
        'X': (lambda x, y: x ^ y, 1),
        'C': (lambda x, y: x & y, 1),
        'D': (lambda x, y: x | y, 1),
        'E': (lambda x, y: 1 - (x ^ y), 1),
        'I': (lambda x, y: (1 - x) | y, 1),
        # Mixing-dimension operators (dims=2), h1[0] vs h2[1]
        "X'": (lambda h1, h2: h1[0] ^ h2[1], 2),
        "C'": (lambda h1, h2: h1[0] & h2[1], 2),
        "D'": (lambda h1, h2: h1[0] | h2[1], 2),
        "E'": (lambda h1, h2: 1 - (h1[0] ^ h2[1]), 2),
        "I'": (lambda h1, h2: (1 - h1[0]) | h2[1], 2),
    }

    unlearnable_operators = []
    
    print("Evaluating learnability of each operator...")
    for name, (func, dims) in operators.items():
        is_learnable = check_learnability(name, func, dims)
        print(f"- Operator '{name}': {'Learnable' if is_learnable else 'Not Learnable'}")
        if not is_learnable:
            unlearnable_operators.append(name)
            
    print("\nBased on the analysis, the operators that cannot be learned are:")
    print(' '.join(sorted(unlearnable_operators)))

if __name__ == '__main__':
    main()