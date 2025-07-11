import numpy as np
from sklearn.linear_model import LogisticRegression

def solve():
    """
    This function analyzes the learnability of various logical relations
    using a specific sentence embedding heuristic and a logistic regression model.
    """
    # Step 1: Generate all possible 2D binary input vectors
    # h1 = [p1, p2], h2 = [q1, q2]
    inputs = []
    for p1 in [0, 1]:
        for p2 in [0, 1]:
            for q1 in [0, 1]:
                for q2 in [0, 1]:
                    h1 = np.array([p1, p2])
                    h2 = np.array([q1, q2])
                    inputs.append((h1, h2))

    # Step 2: Compute heuristic features for all inputs
    # The heuristic is [h1, h2, |h1 - h2|, h1 ⊙ h2]
    X_data = []
    for h1, h2 in inputs:
        h_diff = np.abs(h1 - h2)
        h_prod = h1 * h2
        feature_vector = np.concatenate([h1, h2, h_diff, h_prod])
        X_data.append(feature_vector)
    X_data = np.array(X_data)
    feature_names = ["p1", "p2", "q1", "q2", "|p1-q1|", "|p2-q2|", "p1*q1", "p2*q2"]

    # Step 3: Define all relation operators as functions
    relations = {
        # Element-wise operators (testing on the first dimension)
        "X (p1 XOR q1)":       lambda h1, h2: np.logical_xor(h1[0], h2[0]),
        "C (p1 AND q1)":       lambda h1, h2: np.logical_and(h1[0], h2[0]),
        "D (p1 OR q1)":        lambda h1, h2: np.logical_or(h1[0], h2[0]),
        "E (p1 <=> q1)":       lambda h1, h2: h1[0] == h2[0],
        "I (p1 -> q1)":        lambda h1, h2: np.logical_or(np.logical_not(h1[0]), h2[0]),

        # Mixing dimension operators
        "X' (p1 XOR q2)":      lambda h1, h2: np.logical_xor(h1[0], h2[1]),
        "C' (p1 AND q2)":      lambda h1, h2: np.logical_and(h1[0], h2[1]),
        "D' (p1 OR q2)":       lambda h1, h2: np.logical_or(h1[0], h2[1]),
        "E' (p1 <=> q2)":      lambda h1, h2: h1[0] == h2[1],
        "I' (p1 -> q2)":       lambda h1, h2: np.logical_or(np.logical_not(h1[0]), h2[1])
    }

    # Step 4: Train a classifier for each relation and evaluate
    unlearnable_relations = []
    print("Analyzing linear separability of logical relations...\n")

    for name, func in relations.items():
        # Compute target labels y for the current relation
        y_data = np.array([func(h1, h2) for h1, h2 in inputs], dtype=int)
        
        # Use Logistic Regression to test for linear separability.
        # High C value encourages finding a perfect separator if one exists.
        model = LogisticRegression(C=1e5, solver='liblinear', tol=1e-6)
        
        # Check if there are at least two classes to classify
        if len(np.unique(y_data)) < 2:
            print(f"Relation: {name}\nSkipping: only one class in target labels.\n")
            continue

        model.fit(X_data, y_data)
        
        # An accuracy of 1.0 means the data is linearly separable and thus "learnable".
        accuracy = model.score(X_data, y_data)
        is_learnable = accuracy == 1.0
        
        if not is_learnable:
            unlearnable_relations.append(name.split(" ")[0])

        print(f"Relation: {name}")
        print(f"Accuracy: {accuracy}")
        print(f"Linearly Separable (Learnable): {'Yes' if is_learnable else 'No'}")
        
        # Step 5: Output the learned equation for transparency
        coef = model.coef_[0]
        intercept = model.intercept_[0]
        
        # Building the equation string with each number
        equation_parts = []
        for i, fname in enumerate(feature_names):
            equation_parts.append(f"({coef[i]:+.2f} * {fname})")
        
        equation = " + ".join(equation_parts) + f" + ({intercept:+.2f})"
        print(f"Learned decision function: sign({equation})\n")
        
    print("-" * 50)
    print("Conclusion:")
    print("The model structure is additive across dimensions and cannot learn relations that mix dimensions.")
    print("\nThe following relations CANNOT be learned:")
    unlearnable_relations.sort()
    print(", ".join(unlearnable_relations))

solve()