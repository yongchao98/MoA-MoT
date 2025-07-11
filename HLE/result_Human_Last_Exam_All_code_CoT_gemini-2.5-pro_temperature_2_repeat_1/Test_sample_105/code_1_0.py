import numpy as np
from scipy.optimize import linprog
from itertools import product

def analyze_learnability():
    """
    Analyzes the learnability of various logical operators given a specific
    feature heuristic and a linear classifier model.
    """
    
    # --- Feature Generation Functions ---

    def get_features_elementwise(p1, p2):
        """
        Computes heuristic features for a single dimension. Bias is included.
        Features: [bias, p1, p2, |p1-p2|, p1*p2]
        """
        # [h1, h2, |h1-h2|, h1*h2] and a bias term
        return np.array([1, p1, p2, abs(p1 - p2), p1 * p2])

    def get_features_mixed(p1, q1, p2, q2):
        """
        Computes heuristic features for two dimensions. Bias is included.
        """
        h1 = np.array([p1, q1])
        h2 = np.array([p2, q2])
        # Concatenate [h1, h2, |h1-h2|, h1*h2]
        concatenated = np.concatenate([h1, h2, np.abs(h1 - h2), h1 * h2])
        # Add bias term at the beginning
        return np.insert(concatenated, 0, 1)

    # --- Linear Separability Checker ---

    def check_learnability(feature_func, target_func, input_space, input_vars):
        """
        Checks if a logical operator is learnable (linearly separable).
        Returns True if learnable, False otherwise.
        """
        A_ub = []
        # Get the number of features from the first possible input
        num_features = len(feature_func(*input_space[0]))
        
        for inputs in input_space:
            input_dict = dict(zip(input_vars, inputs))
            target_val = target_func(**input_dict)
            # Convert target {0, 1} to class label {-1, 1}
            y = 1 if target_val else -1
            # Get the feature vector x for the current input
            x = feature_func(*inputs)
            # Add the inequality constraint -y * x @ w <= -1
            A_ub.append(-y * x)

        A_ub = np.array(A_ub)
        # We require y*(w@x) >= 1, which means -y*(w@x) <= -1.
        b_ub = -np.ones(A_ub.shape[0])
        # We only care about feasibility, so the objective function is all zeros.
        c = np.zeros(num_features)

        # Use a linear programming solver to find if a feasible solution exists.
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(None, None))
        
        # A status of 2 indicates the problem is infeasible.
        return res.status != 2

    # --- Operator and Input Space Definitions ---
    
    operators = {
        # Element-wise operators
        'X': (lambda p1, p2, **k: p1 != p2, get_features_elementwise, list(product([0, 1], repeat=2)), ['p1', 'p2']),
        'C': (lambda p1, p2, **k: p1 and p2, get_features_elementwise, list(product([0, 1], repeat=2)), ['p1', 'p2']),
        'D': (lambda p1, p2, **k: p1 or p2, get_features_elementwise, list(product([0, 1], repeat=2)), ['p1', 'p2']),
        'E': (lambda p1, p2, **k: p1 == p2, get_features_elementwise, list(product([0, 1], repeat=2)), ['p1', 'p2']),
        'I': (lambda p1, p2, **k: (not p1) or p2, get_features_elementwise, list(product([0, 1], repeat=2)), ['p1', 'p2']),
        # Mixing dimensions operators
        "X'": (lambda p1, q2, **k: p1 != q2, get_features_mixed, list(product([0, 1], repeat=4)), ['p1', 'q1', 'p2', 'q2']),
        "C'": (lambda p1, q2, **k: p1 and q2, get_features_mixed, list(product([0, 1], repeat=4)), ['p1', 'q1', 'p2', 'q2']),
        "D'": (lambda p1, q2, **k: p1 or q2, get_features_mixed, list(product([0, 1], repeat=4)), ['p1', 'q1', 'p2', 'q2']),
        "E'": (lambda p1, q2, **k: p1 == q2, get_features_mixed, list(product([0, 1], repeat=4)), ['p1', 'q1', 'p2', 'q2']),
        "I'": (lambda p1, q2, **k: (not p1) or q2, get_features_mixed, list(product([0, 1], repeat=4)), ['p1', 'q1', 'p2', 'q2']),
    }
    
    unlearnable_operators = []
    
    print("--- Analysis Results ---")
    for name, (func, feat_func, space, in_vars) in operators.items():
        is_learnable = check_learnability(feat_func, func, space, in_vars)
        print(f"Operator '{name}': Learnable = {is_learnable}")
        if not is_learnable:
            unlearnable_operators.append(name)
            
    print("\n-------------------------------------------")
    print("Conclusion: The operators that can not be learned are:")
    for op in unlearnable_operators:
        print(f"- {op}")
    print("\nCombining these gives the string:", "".join(sorted(unlearnable_operators)))

if __name__ == '__main__':
    analyze_learnability()