import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product
import warnings

# Suppress convergence warnings for this specific case
warnings.filterwarnings('ignore', category=UserWarning)

def heuristic(h1, h2):
    """Computes the heuristic [h1, h2, |h1-h2|, h1*h2]"""
    h1 = np.array(h1)
    h2 = np.array(h2)
    diff = np.abs(h1 - h2)
    prod = h1 * h2
    return np.concatenate([h1, h2, diff, prod])

# --- Define Logical Operators ---
# Note: h1 and h2 are tuples representing embeddings, e.g., (0,) or (1, 0)
def op_X(h1, h2): return h1[0] ^ h2[0]
def op_C(h1, h2): return h1[0] & h2[0]
def op_D(h1, h2): return h1[0] | h2[0]
def op_E(h1, h2): return 1 if h1[0] == h2[0] else 0
def op_I(h1, h2): return 1 if not h1[0] or h2[0] else 0

# Mixing-dimension operators use h1[0] and h2[1]
def op_X_prime(h1, h2): return h1[0] ^ h2[1]
def op_C_prime(h1, h2): return h1[0] & h2[1]
def op_D_prime(h1, h2): return h1[0] | h2[1]
def op_E_prime(h1, h2): return 1 if h1[0] == h2[1] else 0
def op_I_prime(h1, h2): return 1 if not h1[0] or h2[1] else 0

def test_operator(op_func, dim):
    """
    Tests if a logical operator is learnable.
    Returns True if learnable (100% accuracy), False otherwise.
    """
    possible_values = list(product([0, 1], repeat=dim))
    all_inputs = list(product(possible_values, repeat=2)) # (h1, h2) tuples

    X_data = np.array([heuristic(h1, h2) for h1, h2 in all_inputs])
    y_labels = np.array([op_func(h1, h2) for h1, h2 in all_inputs])
    
    # If all labels are the same, it's trivially learnable
    if len(np.unique(y_labels)) < 2:
        return True

    # Use a large C for minimal regularization to see if a separating plane exists.
    # liblinear is a good solver for small datasets.
    model = LogisticRegression(random_state=0, solver='liblinear', C=1e5)
    model.fit(X_data, y_labels)
    
    return model.score(X_data, y_labels) == 1.0

def find_unlearnable_operators():
    """
    Main function to find and list all unlearnable operators.
    """
    operators_to_test = {
        'X': (op_X, 1), 'C': (op_C, 1), 'D': (op_D, 1), 'E': (op_E, 1), 'I': (op_I, 1),
        "X'": (op_X_prime, 2), "C'": (op_C_prime, 2), "D'": (op_D_prime, 2),
        "E'": (op_E_prime, 2), "I'": (op_I_prime, 2)
    }

    not_learnable = []
    print("Running analysis...")
    for name, (func, dim) in operators_to_test.items():
        if not test_operator(func, dim):
            not_learnable.append(name)
            
    # Sort for consistent output, e.g., X'E' not E'X'
    not_learnable.sort()
    
    print("\nBased on the analysis, the operators that cannot be learned are:")
    print("".join(not_learnable))

if __name__ == '__main__':
    find_unlearnable_operators()
