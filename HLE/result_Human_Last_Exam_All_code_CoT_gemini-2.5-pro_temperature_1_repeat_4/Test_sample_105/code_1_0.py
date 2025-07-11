import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product
import warnings

# Suppress convergence warnings for cleaner output
warnings.filterwarnings("ignore")

def solve():
    """
    Determines which logical operators are not learnable by a logistic regression
    model using the given heuristic representation.
    """
    unlearnable_ops = []
    learnable_ops = []
    
    # 1. Test Element-wise Operators (X, C, D, E, I)
    # These operate on a single dimension. We have 4 possible input pairs (0,0), (0,1), (1,0), (1,1).
    x = np.array([0, 0, 1, 1]) # Represents h1[i] for the 4 cases
    y = np.array([0, 1, 0, 1]) # Represents h2[i] for the 4 cases

    # Heuristic features: [h1, h2, |h1-h2|, h1*h2]
    H_elementwise = np.vstack([x, y, np.abs(x - y), x * y]).T

    # Define the target outputs for each element-wise operator
    targets_elementwise = {
        'X (XOR)': (x != y).astype(int),
        'C (Conjunction)': (x * y),
        'D (Disjunction)': np.minimum(1, x + y),
        'E (Equivalence)': (x == y).astype(int),
        'I (Implication)': ( (1 - x) | y ), # not x or y
    }

    print("--- Analyzing Element-wise Operators ---")
    for name, y_target in targets_elementwise.items():
        if len(np.unique(y_target)) > 1: # Check if the problem is non-trivial
            model = LogisticRegression(solver='liblinear', C=1e9) # Use high C to minimize regularization
            model.fit(H_elementwise, y_target)
            score = model.score(H_elementwise, y_target)
            if score < 1.0:
                unlearnable_ops.append(name.split(' ')[0])
            else:
                learnable_ops.append(name.split(' ')[0])
        else: # Trivial case
            learnable_ops.append(name.split(' ')[0])

    # 2. Test Mixing-Dimensions Operators (X', C', D', E', I')
    # These operate on two dimensions, h1[i] and h2[j].
    # We need all 2^4 = 16 combinations for (h1[i], h2[i], h1[j], h2[j])
    # Let's use x1, y1 for dim i and x2, y2 for dim j.
    inputs = list(product([0, 1], repeat=4)) # Generates (x1, y1, x2, y2) tuples
    
    H_mix = []
    # The target depends on h1[i] (x1) and h2[j] (y2)
    y_mix_targets = {
        "X' (XOR mixing)": [],
        "C' (Conjunction mixing)": [],
        "D' (Disjunction mixing)": [],
        "E' (Equivalence mixing)": [],
        "I' (Implication mixing)": [],
    }

    for x1, y1, x2, y2 in inputs:
        # Heuristic features for dimension i and j are concatenated
        feature_vec = [x1, y1, abs(x1 - y1), x1 * y1, x2, y2, abs(x2 - y2), x2 * y2]
        H_mix.append(feature_vec)
        
        # Define target outputs which mix dimensions (x1 and y2)
        y_mix_targets["X' (XOR mixing)"].append(int(x1 != y2))
        y_mix_targets["C' (Conjunction mixing)"].append(x1 * y2)
        y_mix_targets["D' (Disjunction mixing)"].append(int(x1 or y2))
        y_mix_targets["E' (Equivalence mixing)"].append(int(x1 == y2))
        y_mix_targets["I' (Implication mixing)"].append(int((1 - x1) or y2))

    H_mix = np.array(H_mix)

    print("\n--- Analyzing Mixing-Dimensions Operators ---")
    for name, y_target in y_mix_targets.items():
        y_target = np.array(y_target)
        if len(np.unique(y_target)) > 1:
            model = LogisticRegression(solver='liblinear', C=1e9)
            model.fit(H_mix, y_target)
            score = model.score(H_mix, y_target)
            if score < 1.0:
                unlearnable_ops.append(name.split(' ')[0])
            else:
                learnable_ops.append(name.split(' ')[0])
        else: # Trivial case
            learnable_ops.append(name.split(' ')[0])
            
    print("\n--- Conclusion ---")
    print(f"Learnable operators: {sorted(learnable_ops)}")
    print(f"Unlearnable operators: {sorted(unlearnable_ops)}")
    
    # Format the final list to match the answer choices
    final_unlearnable_list = "".join(sorted([op for op in unlearnable_ops if "'" in op]))
    print(f"\nThe list of unlearnable operators from the choices is: {final_unlearnable_list}")


solve()
<<<I>>>