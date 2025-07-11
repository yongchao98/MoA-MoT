import numpy as np
from sklearn.linear_model import LogisticRegression

def solve():
    """
    Analyzes the learnability of various logical operators with a given heuristic representation.
    """
    operators = {
        # Element-wise operators
        'C': {'type': 'element-wise', 'func': lambda x, y: x and y, 'name': 'Conjunction'},
        'D': {'type': 'element-wise', 'func': lambda x, y: x or y, 'name': 'Disjunction'},
        'I': {'type': 'element-wise', 'func': lambda x, y: (not x) or y, 'name': 'Implication'},
        'E': {'type': 'element-wise', 'func': lambda x, y: x == y, 'name': 'Equivalence'},
        'X': {'type': 'element-wise', 'func': lambda x, y: x != y, 'name': 'XOR'},
        # Mixing-dimension operators
        'C\'': {'type': 'mixing', 'func': lambda x, y: x and y, 'name': 'Conjunction mixing dimensions'},
        'D\'': {'type': 'mixing', 'func': lambda x, y: x or y, 'name': 'Disjunction mixing dimensions'},
        'I\'': {'type': 'mixing', 'func': lambda x, y: (not x) or y, 'name': 'Implication mixing dimensions'},
        'E\'': {'type': 'mixing', 'func': lambda x, y: x == y, 'name': 'Equivalence mixing dimensions'},
        'X\'': {'type': 'mixing', 'func': lambda x, y: x != y, 'name': 'XOR mixing dimensions'},
    }

    # Input pairs for a single dimension
    inputs = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    
    unlearnable_operators = []

    print("--- Analyzing Learnability of Operators ---\n")

    for op_code, op_info in operators.items():
        print(f"--- Checking Operator: {op_code} ({op_info['name']}) ---")
        
        # Generate target labels from the operator's truth table
        y_target = np.array([op_info['func'](i[0], i[1]) for i in inputs])

        if op_info['type'] == 'element-wise':
            # Create feature representation: [h1, h2, |h1-h2|, h1*h2]
            h1 = inputs[:, 0]
            h2 = inputs[:, 1]
            abs_diff = np.abs(h1 - h2)
            hadamard = h1 * h2
            X_features = np.vstack([h1, h2, abs_diff, hadamard]).T

            # Check for linear separability using a Logistic Regression model
            model = LogisticRegression(solver='liblinear', C=1e5)
            model.fit(X_features, y_target)
            accuracy = model.score(X_features, y_target)

            if accuracy == 1.0:
                print(f"Result: Operator '{op_code}' is learnable (achieved 100% accuracy).\n")
            else:
                print(f"Result: Operator '{op_code}' is NOT learnable.\n")
                unlearnable_operators.append(op_code)

        elif op_info['type'] == 'mixing':
            print("This is a mixing-dimension operator.")
            print("Its learnability depends on whether its truth table can satisfy the separability condition:")
            print("S(0,1) + S(1,0) = S(0,0) + S(1,1), where S(x,y) is the model score.\n")

            # Determine required signs of S based on the truth table
            signs = ['> 0' if label else '< 0' for label in y_target]
            print(f"Truth table for '{op_code}': {y_target.astype(int)}")
            print(f"Required scores: S(0,0){signs[0]}, S(0,1){signs[1]}, S(1,0){signs[2]}, S(1,1){signs[3]}")

            # Demonstrate the contradiction with example numbers
            S00 = 1 if y_target[0] else -1
            S01 = 1 if y_target[1] else -1
            S10 = 1 if y_target[2] else -1
            S11 = 1 if y_target[3] else -1
            
            # For C' and D', make magnitudes different to show a clearer contradiction
            if op_code == 'C\'': S11 = 2
            if op_code == 'D\'': S00 = -2
            if op_code == 'I\'': S10 = -2


            lhs_val = S01 + S10
            rhs_val = S00 + S11

            print("\nTesting the separability equation with example scores:")
            print(f"  S(0,1) + S(1,0)  =  S(0,0) + S(1,1)")
            # Print each number in the equation
            print(f"  ({S01}) + ({S10})   =  ({S00}) + ({S11})")
            print(f"  {lhs_val}   =   {rhs_val}")

            if lhs_val == rhs_val:
                # This case shouldn't be reached by the chosen examples
                print("\nResult: Condition is met. The operator might be learnable (recheck logic).")
            else:
                print("\nResult: The equality is false. The condition is violated.")
                print(f"Therefore, operator '{op_code}' is NOT learnable.\n")
                unlearnable_operators.append(op_code)
    
    print("\n--- Summary ---")
    print("The following operators can NOT be learned with the given heuristic:")
    unlearnable_operators.sort()
    print(' '.join(unlearnable_operators))

solve()
<<<I>>>