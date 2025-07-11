import numpy as np
from itertools import product

def get_feature_vector(h1, h2):
    """
    Computes the 4-dimensional feature vector based on the weight-sharing model.
    """
    h1 = np.array(h1)
    h2 = np.array(h2)
    
    # Ensure inputs are numpy arrays for vectorized operations
    sum_h1 = np.sum(h1)
    sum_h2 = np.sum(h2)
    sum_abs_diff = np.sum(np.abs(h1 - h2))
    sum_hadamard = np.sum(h1 * h2)
    
    return (sum_h1, sum_h2, sum_abs_diff, sum_hadamard)

def check_learnability(operator_name, target_function, dims):
    """
    Checks if an operator is learnable by looking for collisions.
    A collision occurs if the same feature vector maps to different target labels.
    """
    feature_map = {}
    
    # Generate all possible input embeddings h1 and h2 of size 'dims'
    possible_vectors = list(product([0, 1], repeat=dims))
    
    for h1 in possible_vectors:
        for h2 in possible_vectors:
            
            # 1. Compute the feature vector for this (h1, h2) pair
            features = get_feature_vector(h1, h2)
            
            # 2. Compute the target label for this (h1, h2) pair
            label = target_function(h1, h2)
            
            # 3. Check for collisions
            if features not in feature_map:
                feature_map[features] = label
            elif feature_map[features] != label:
                # Collision detected! The same features require different labels.
                return False # Not learnable
                
    return True # Learnable

def main():
    # --- Define logical operators ---

    # Element-wise operators (use dims=1)
    # The 'h' vectors will have one element, e.g., h1 = (p1,) h2 = (p2,)
    elementwise_ops = {
        "C": lambda h1, h2: h1[0] & h2[0],
        "D": lambda h1, h2: h1[0] | h2[0],
        "X": lambda h1, h2: h1[0] ^ h2[0],
        "E": lambda h1, h2: 1 if h1[0] == h2[0] else 0,
        "I": lambda h1, h2: 1 if not h1[0] or h2[0] else 0,
    }

    # Mixed-dimension operators (use dims=2)
    # The 'h' vectors will have two elements, e.g., h1=(p1, q1), h2=(p2, q2)
    # We define the relation between the first dim of h1 (p1) and second dim of h2 (q2)
    mixed_ops = {
        "C'": lambda h1, h2: h1[0] & h2[1],
        "D'": lambda h1, h2: h1[0] | h2[1],
        "X'": lambda h1, h2: h1[0] ^ h2[1],
        "E'": lambda h1, h2: 1 if h1[0] == h2[1] else 0,
        "I'": lambda h1, h2: 1 if not h1[0] or h2[1] else 0,
    }
    
    unlearnable = []
    
    print("Analyzing Element-wise Operators (dims=1)...")
    for name, func in elementwise_ops.items():
        if not check_learnability(name, func, dims=1):
            unlearnable.append(name)

    print("Analyzing Mixed-Dimension Operators (dims=2)...")
    for name, func in mixed_ops.items():
        if not check_learnability(name, func, dims=2):
            unlearnable.append(name)
            
    # The problem asks to print the numbers/results in an equation format
    print("\nResult:")
    print("Unlearnable Operators = { " + ", ".join(sorted(unlearnable)) + " }")

if __name__ == "__main__":
    main()