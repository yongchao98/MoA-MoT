import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
import warnings

# Suppress convergence warnings for this demonstration
warnings.filterwarnings("ignore", category=ConvergenceWarning)

def solve():
    """
    This function programmatically determines which logical operators cannot be learned
    by a logistic regression model on the given heuristic feature representation.
    """
    # 1. Generate all 16 possible inputs for 2D binary embeddings
    inputs = []
    for p1 in [0, 1]:
        for q1 in [0, 1]:
            for p2 in [0, 1]:
                for q2 in [0, 1]:
                    inputs.append({'p1': p1, 'q1': q1, 'p2': p2, 'q2': q2})

    # 2. For each input, compute the heuristic feature vector
    # The heuristic is [h1, h2, |h1-h2|, h1*h2]
    X_features = []
    for i in inputs:
        h1 = np.array([i['p1'], i['q1']])
        h2 = np.array([i['p2'], i['q2']])
        abs_diff = np.abs(h1 - h2)
        h_prod = h1 * h2
        feature_vector = np.concatenate([h1, h2, abs_diff, h_prod])
        X_features.append(feature_vector)
    X = np.array(X_features)

    # 3. Define the target logical operators
    # Note: For element-wise operators, we test on the first dimension (p1, p2) without loss of generality.
    operators = {
        'X': lambda i: i['p1'] ^ i['p2'],
        'C': lambda i: i['p1'] & i['p2'],
        'D': lambda i: i['p1'] | i['p2'],
        'E': lambda i: int(i['p1'] == i['p2']),
        'I': lambda i: int(not i['p1'] or i['p2']),
        'X\'': lambda i: i['p1'] ^ i['q2'],
        'C\'': lambda i: i['p1'] & i['q2'],
        'D\'': lambda i: i['p1'] | i['q2'],
        'E\'': lambda i: int(i['p1'] == i['q2']),
        'I\'': lambda i: int(not i['p1'] or i['q2']),
    }

    not_learnable = []
    print("Evaluating operators for learnability...")
    print("-" * 40)

    # 4. For each operator, train a model and check if it's perfectly learnable
    for name, op_func in operators.items():
        y = np.array([op_func(i) for i in inputs])

        # A problem is only interesting if there are two classes to separate
        if len(np.unique(y)) < 2:
            print(f"Operator '{name}': Learnable (trivial case, single class).")
            continue
        
        # Use a high C value to find a separating hyperplane if one exists
        model = LogisticRegression(solver='liblinear', C=1e5, tol=1e-8, max_iter=1000)
        model.fit(X, y)
        accuracy = model.score(X, y)

        if accuracy < 1.0:
            not_learnable.append(name)
            print(f"Operator '{name}': NOT Learnable (Accuracy: {accuracy:.2f})")
        else:
            print(f"Operator '{name}': Learnable (Accuracy: {accuracy:.2f})")
            
    # 5. Output the final list
    print("-" * 40)
    print("The following operators cannot be learned with the given heuristic:")
    # The strange request in the prompt: "output each number in the final equation!"
    # is interpreted as printing the final list of non-learnable operators.
    result_string = "".join(sorted(not_learnable))
    print(result_string)
    
solve()
<<<I>>>