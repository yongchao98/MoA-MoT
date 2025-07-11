import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product

def solve():
    """
    Determines which logical relations are not learnable by a logistic regression
    model using the specified feature heuristic.
    """

    # 1. Define the heuristic function
    def heuristic(h1, h2):
        """
        Computes the heuristic [h1, h2, |h1-h2|, h1*h2]
        """
        h1 = np.array(h1)
        h2 = np.array(h2)
        diff = np.abs(h1 - h2)
        hadamard = h1 * h2
        return np.concatenate([h1, h2, diff, hadamard])

    # 2. Define the logical relations
    relations = {
        # Element-wise relations (use 1 dimension)
        'X': {'func': lambda h1, h2: h1[0] ^ h2[0], 'dims': 1},
        'C': {'func': lambda h1, h2: h1[0] & h2[0], 'dims': 1},
        'D': {'func': lambda h1, h2: h1[0] | h2[0], 'dims': 1},
        'E': {'func': lambda h1, h2: int(h1[0] == h2[0]), 'dims': 1},
        'I': {'func': lambda h1, h2: int(not h1[0] or h2[0]), 'dims': 1},
        
        # Mixing-dimension relations (use 2 dimensions)
        # The relation is between the first element of h1 and the second of h2
        'X\'': {'func': lambda h1, h2: h1[0] ^ h2[1], 'dims': 2},
        'C\'': {'func': lambda h1, h2: h1[0] & h2[1], 'dims': 2},
        'D\'': {'func': lambda h1, h2: h1[0] | h2[1], 'dims': 2},
        'E\'': {'func': lambda h1, h2: int(h1[0] == h2[1]), 'dims': 2},
        'I\'': {'func': lambda h1, h2: int(not h1[0] or h2[1]), 'dims': 2},
    }

    unlearnable_relations = []
    
    # Order to match choice 'I'
    relation_order = ['X', 'C', 'D', 'E', 'I', 'X\'', 'C\'', 'D\'', 'E\'', 'I\'']

    # 3-6. Generate data, transform features, train, and evaluate for each relation
    for name in relation_order:
        props = relations[name]
        dims = props['dims']
        relation_func = props['func']

        possible_vectors = list(product([0, 1], repeat=dims))
        
        X_data = []
        y_data = []

        for h1 in possible_vectors:
            for h2 in possible_vectors:
                features = heuristic(h1, h2)
                X_data.append(features)
                label = relation_func(h1, h2)
                y_data.append(label)

        X = np.array(X_data)
        y = np.array(y_data)

        if len(np.unique(y)) < 2:
            learnable = True
        else:
            # Use a high C value to disable regularization, focusing on separability
            model = LogisticRegression(solver='liblinear', C=1e9, random_state=0)
            model.fit(X, y)
            score = model.score(X, y)
            learnable = (score == 1.0)
        
        if not learnable:
            unlearnable_relations.append(name)

    # 7. Collect and Report
    final_answer_string = "".join(unlearnable_relations)
    print("List of relation operators that cannot be learned:")
    print(final_answer_string)

solve()