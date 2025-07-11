import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product

def get_heuristic_features(h1, h2):
    """
    Computes the heuristic [h1, h2, |h1-h2|, h1*h2]
    """
    if h1.ndim == 1:
        h1 = h1.reshape(1, -1)
        h2 = h2.reshape(1, -1)
        
    abs_diff = np.abs(h1 - h2)
    hadamard_prod = h1 * h2
    # Concatenate all features for each sample
    features = np.concatenate([h1, h2, abs_diff, hadamard_prod], axis=1)
    return features

# We need at least 2 dimensions to test mixing-dimension operators
DIM = 2

# Generate all 2^(2*DIM) = 16 possible inputs
# h1 = (p1, q1), h2 = (p2, q2)
possible_values = [0, 1]
all_inputs = list(product(possible_values, repeat=2*DIM))

# Reshape into (h1, h2) pairs
X_source = []
for p in all_inputs:
    h1 = np.array(p[:DIM])
    h2 = np.array(p[DIM:])
    X_source.append((h1, h2))

# Generate heuristic features for all inputs
X_features = np.array([get_heuristic_features(h1, h2)[0] for h1, h2 in X_source])

# Define all logical operators to test
operators = {
    # Element-wise operators (test on first dimension)
    'X': lambda h1, h2: h1[0] ^ h2[0],
    'C': lambda h1, h2: h1[0] & h2[0],
    'D': lambda h1, h2: h1[0] | h2[0],
    'E': lambda h1, h2: int(h1[0] == h2[0]),
    'I': lambda h1, h2: int(not h1[0] or h2[0]),
    
    # Mixing-dimension operators (use dim 0 and 1)
    'X\'': lambda h1, h2: h1[0] ^ h2[1],
    'C\'': lambda h1, h2: h1[0] & h2[1],
    'D\'': lambda h1, h2: h1[0] | h2[1],
    'E\'': lambda h1, h2: int(h1[0] == h2[1]),
    'I\'': lambda h1, h2: int(not h1[0] or h2[1])
}

unlearnable_operators = []

print("Analyzing operators for learnability...\n")

for name, op_func in operators.items():
    # Generate target labels y for the current operator
    y_target = np.array([op_func(h1, h2) for h1, h2 in X_source])
    
    # If all labels are the same, it's trivially learnable.
    # We are interested in non-trivial cases.
    if len(np.unique(y_target)) < 2:
        print(f"Operator '{name}': Trivially learnable (constant output).")
        continue

    # Train a logistic regression model.
    # High C means low regularization, to fit the data if possible.
    model = LogisticRegression(C=1e9, solver='liblinear', tol=1e-6)
    model.fit(X_features, y_target)
    
    # Check if the model learned the function perfectly
    accuracy = model.score(X_features, y_target)
    
    print(f"Operator '{name}': Accuracy = {accuracy:.2f}. {'Learnable' if accuracy == 1.0 else 'Not Learnable'}")
    if accuracy < 1.0:
        unlearnable_operators.append(name)

print("\n-------------------------------------------")
print(f"Conclusion: The set of relations that can not be learned is:")
print(unlearnable_operators)
print("-------------------------------------------")
<<<I>>>