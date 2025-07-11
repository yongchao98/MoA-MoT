import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product

def get_heuristic_features(h1, h2):
    """
    Computes the heuristic representation for one or more dimensions.
    h1 and h2 are numpy arrays of the same length.
    """
    abs_diff = np.abs(h1 - h2)
    hadamard = h1 * h2
    return np.concatenate([h1, h2, abs_diff, hadamard])

def check_learnability(target_func, is_mixing):
    """
    Checks if an operator is learnable by the model.
    An operator is learnable if a logistic regression model can achieve 100% 
    accuracy on all possible inputs, meaning the classes are linearly separable
    in the heuristic feature space.
    """
    if not is_mixing:
        # Element-wise operator, check on a single dimension.
        # The inputs are all possible pairs of (h1[i], h2[i]).
        inputs = np.array(list(product([0, 1], repeat=2)))
        h1_vals = inputs[:, 0]
        h2_vals = inputs[:, 1]
        
        # Generate feature vectors X for each input pair.
        # Each feature vector has 4 elements.
        X = np.array([get_heuristic_features(np.array([h1]), np.array([h2])) for h1, h2 in zip(h1_vals, h2_vals)])
        
        # Generate target labels y.
        y = target_func(h1_vals, h2_vals)
        
    else:
        # Mixing-dimension operator, check on two dimensions (i, j).
        # The inputs are all 16 possible tuples of (h1[i], h1[j], h2[i], h2[j]).
        inputs = np.array(list(product([0, 1], repeat=4)))
        h1i = inputs[:, 0]
        h1j = inputs[:, 1]
        h2i = inputs[:, 2]
        h2j = inputs[:, 3]
        
        # Generate feature vectors X. For each of the 16 inputs, h1 is (h1i, h1j)
        # and h2 is (h2i, h2j). The resulting feature vector has 8 elements.
        X_list = []
        for i in range(inputs.shape[0]):
            h1 = np.array([h1i[i], h1j[i]])
            h2 = np.array([h2i[i], h2j[i]])
            X_list.append(get_heuristic_features(h1, h2))
        X = np.array(X_list)
        
        # Generate target labels y. The target function for mixing operators
        # is defined on h1[i] and h2[j].
        y = target_func(h1i, h2j)

    # If all labels are the same, the function is constant and trivially learnable.
    if len(np.unique(y)) < 2:
        return True
        
    # Train a logistic regression model.
    # We use a large C to minimize regularization, forcing the model to find a
    # separating hyperplane if one exists.
    model = LogisticRegression(solver='liblinear', C=1e9, random_state=0)
    model.fit(X, y)
    
    # Check if the model learned the training data perfectly.
    accuracy = model.score(X, y)
    
    return accuracy == 1.0

def main():
    """
    Main function to define operators, test their learnability, and print the result.
    """
    # Define all operators to be tested.
    # Each entry is a tuple: (target_function, is_mixing_dimension_operator)
    operators = {
        'X (element-wise XOR)': (lambda h1, h2: np.logical_xor(h1, h2), False),
        'C (element-wise Conjunction)': (lambda h1, h2: np.logical_and(h1, h2), False),
        'D (element-wise Disjunction)': (lambda h1, h2: np.logical_or(h1, h2), False),
        'E (element-wise Equivalence)': (lambda h1, h2: np.equal(h1, h2), False),
        'I (element-wise Implication)': (lambda h1, h2: np.logical_or(np.logical_not(h1), h2), False),
        "X' (XOR mixing dimensions)": (lambda h1i, h2j: np.logical_xor(h1i, h2j), True),
        "C' (Conjunction mixing dimensions)": (lambda h1i, h2j: np.logical_and(h1i, h2j), True),
        "D' (Disjunction mixing dimensions)": (lambda h1i, h2j: np.logical_or(h1i, h2j), True),
        "E' (Equivalence mixing dimensions)": (lambda h1i, h2j: np.equal(h1i, h2j), True),
        "I' (Implication mixing dimensions)": (lambda h1i, h2j: np.logical_or(np.logical_not(h1i), h2j), True),
    }

    print("Analyzing learnability of logical operators...")
    not_learnable_operators = []
    
    for name, (func, is_mixing) in operators.items():
        if not check_learnability(func, is_mixing):
            # Extract the short name (e.g., "X'") from the full name
            short_name = name.split(' ')[0]
            not_learnable_operators.append(short_name)

    print("\nBased on the analysis, the relation operators that CANNOT be learned are:")
    # Sort the list for a consistent, clean output that matches the answer choices.
    print(' '.join(sorted(not_learnable_operators)))

if __name__ == "__main__":
    main()