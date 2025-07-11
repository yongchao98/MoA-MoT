import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product

def solve():
    """
    This function tests the learnability of various logical relations with a specific
    heuristic sentence embedding composition method.
    
    The method is learnable if a logistic regression model can achieve 100% accuracy 
    on all possible inputs.
    """
    
    # Define the logical relations to be tested.
    # For element-wise relations, we check if the relation holds for ANY dimension (OR aggregation).
    # h1 = [p1, q1], h2 = [p2, q2]
    relations = {
        'X': lambda h1, h2: (h1[0] ^ h2[0]) or (h1[1] ^ h2[1]),
        'C': lambda h1, h2: (h1[0] & h2[0]) or (h1[1] & h2[1]),
        'D': lambda h1, h2: (h1[0] | h2[0]) or (h1[1] | h2[1]),
        'E': lambda h1, h2: (h1[0] == h2[0]) or (h1[1] == h2[1]),
        'I': lambda h1, h2: ((1 - h1[0]) | h2[0]) or ((1 - h1[1]) | h2[1]),
        'X\'': lambda h1, h2: h1[0] ^ h2[1],
        'C\'': lambda h1, h2: h1[0] & h2[1],
        'D\'': lambda h1, h2: h1[0] | h2[1],
        'E\'': lambda h1, h2: h1[0] == h2[1],
        'I\'': lambda h1, h2: (1 - h1[0]) | h2[1],
    }

    unlearnable_relations = []

    # Generate all possible 2-dim binary vector inputs
    possible_values = [0, 1]
    all_inputs = list(product(possible_values, repeat=4))

    for name, func in relations.items():
        X_data = []
        y_data = []

        for p1, q1, p2, q2 in all_inputs:
            h1 = np.array([p1, q1])
            h2 = np.array([p2, q2])

            # Compute the heuristic: [h1, h2, |h1-h2|, h1*h2]
            h_abs_diff = np.abs(h1 - h2)
            h_prod = h1 * h2
            
            feature_vector = np.concatenate([h1, h2, h_abs_diff, h_prod])
            X_data.append(feature_vector)
            
            # Compute the target label
            y_data.append(int(func(h1, h2)))
        
        X_data = np.array(X_data)
        y_data = np.array(y_data)

        # A relation is trivially learnable if it's always true or always false.
        # We only need to test non-trivial cases.
        if len(np.unique(y_data)) < 2:
            continue
            
        # Use a logistic regression model to test for linear separability.
        # A high C value minimizes regularization to focus on finding a separating hyperplane.
        model = LogisticRegression(solver='liblinear', C=1e5, random_state=42)
        model.fit(X_data, y_data)
        
        # If the score is not 1.0, the data is not linearly separable by the model.
        score = model.score(X_data, y_data)
        
        if score < 1.0:
            unlearnable_relations.append(name)
            
    print("The following relations cannot be learned by the model:")
    # Sort for consistent ordering, matching option I
    print(sorted(unlearnable_relations, key=lambda x: "XCDEI".find(x[0])))


solve()
<<<I>>>