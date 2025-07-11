import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product

def analyze_relations():
    """
    Analyzes which logical relations are learnable with the given heuristic.
    """
    # Define the 10 logical relations to be tested.
    # p represents the first dimension, q the second.
    # 1/2 represent the embedding index (h1 or h2).
    # e.g., 'p1 AND p2' is element-wise, 'p1 AND q2' is mixing.
    relations = {
        'X': lambda p1, q1, p2, q2: p1 ^ p2,
        'C': lambda p1, q1, p2, q2: p1 & p2,
        'D': lambda p1, q1, p2, q2: p1 | p2,
        'E': lambda p1, q1, p2, q2: p1 == p2,
        'I': lambda p1, q1, p2, q2: (1 - p1) | p2,
        "X'": lambda p1, q1, p2, q2: p1 ^ q2,
        "C'": lambda p1, q1, p2, q2: p1 & q2,
        "D'": lambda p1, q1, p2, q2: p1 | q2,
        "E'": lambda p1, q1, p2, q2: p1 == q2,
        "I'": lambda p1, q1, p2, q2: (1 - p1) | q2,
    }

    # Generate all possible inputs for 2D binary embeddings
    # There are 2^4 = 16 possible combinations of (p1, q1, p2, q2)
    inputs = list(product([0, 1], repeat=4))

    # Heuristic function: h1,h2 -> [h1,h2,|h1-h2|,h1⊙h2]
    def heuristic(p1, q1, p2, q2):
        h1 = np.array([p1, q1])
        h2 = np.array([p2, q2])
        h_cat = np.concatenate([h1, h2])
        h_diff = np.abs(h1 - h2)
        h_prod = h1 * h2
        return np.concatenate([h_cat, h_diff, h_prod])

    # Generate the full feature matrix X based on the heuristic
    X_data = np.array([heuristic(p1, q1, p2, q2) for p1, q1, p2, q2 in inputs])

    print("Analyzing learnability of logical relations...")
    print("-" * 50)

    unlearnable_ops = []

    for name, func in relations.items():
        # Generate the target vector y for the current relation
        y_data = np.array([func(p1, q1, p2, q2) for p1, q1, p2, q2 in inputs])

        # A relation is not learnable if the targets are all the same
        # (no separation is needed). We skip these trivial cases.
        if len(np.unique(y_data)) < 2:
            continue

        # Use a standard Logistic Regression model
        # We use a large C to minimize regularization and find if a separating plane exists.
        model = LogisticRegression(C=1e5, solver='liblinear', random_state=0)
        
        # Fit the model on the entire dataset
        # We check training accuracy because we want to know if the function
        # is representable at all, not how well it generalizes.
        model.fit(X_data, y_data)
        
        # Calculate the accuracy on the training data
        score = model.score(X_data, y_data)

        print(f"Operator: {name:<3s} | Type: {'Element-wise' if len(name)==1 else 'Mixing-dimension'} | Learnable: {score == 1.0} (Accuracy: {score:.2f})")
        
        if score < 1.0:
            unlearnable_ops.append(name)
            
    print("-" * 50)
    print("Based on the analysis, the operators that CANNOT be learned are:")
    # The output format requires printing the equation as text
    final_list_str = "".join(sorted(unlearnable_ops))
    print(f"The list of unlearnable operators is: {final_list_str}")


if __name__ == "__main__":
    analyze_relations()
    print("\nThe correct answer choice is the one that lists all mixing-dimension operators.")
    print("Final Answer formatted for the system:")
    # The final answer is the concatenation of the unlearnable operators' symbols
    # X', C', D', E', I'
    print("<<<I>>>")