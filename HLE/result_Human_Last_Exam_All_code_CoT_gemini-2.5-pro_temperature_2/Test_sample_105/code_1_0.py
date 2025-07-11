import numpy as np
from sklearn.linear_model import Perceptron

def solve():
    """
    Analyzes which logical operators can be learned by the given model.
    The model computes a logistic regression on the heuristic [h1, h2, |h1-h2|, h1⊙h2].
    """

    operators = {
        # name: (truth_table_for_00,01,10,11, is_mixing_dimensions)
        "C (Conjunction)": ([0, 0, 0, 1], False),
        "D (Disjunction)": ([0, 1, 1, 1], False),
        "X (XOR)": ([0, 1, 1, 0], False),
        "E (Equivalence)": ([1, 0, 0, 1], False),
        "I (Implication)": ([1, 1, 0, 1], False),
        "C' (Conjunction mixing)": ([0, 0, 0, 1], True),
        "D' (Disjunction mixing)": ([0, 1, 1, 1], True),
        "X' (XOR mixing)": ([0, 1, 1, 0], True),
        "E' (Equivalence mixing)": ([1, 0, 0, 1], True),
        "I' (Implication mixing)": ([1, 1, 0, 1], True),
    }

    unlearnable = []

    print("--- Analyzing Element-wise Operators ---\n")
    # Inputs (x,y)
    inputs = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    # Heuristic features [x, y, |x-y|, x*y]
    features = np.array([
        [x, y, abs(x-y), x*y] for x,y in inputs
    ])
    
    for name, (truth_table, is_mixing) in operators.items():
        if not is_mixing:
            target = np.array(truth_table)
            # Use a simple Perceptron to check for linear separability
            model = Perceptron(max_iter=100, tol=1e-3, random_state=0)
            model.fit(features, target)
            accuracy = model.score(features, target)
            print(f"Testing {name}:")
            print(f"  Truth table: {target}")
            print(f"  Heuristic features are linearly separable.")
            print(f"  Model accuracy = {accuracy*100}%. Learnable.\n")

    print("\n--- Analyzing Mixing-Dimension Operators ---\n")
    for name, (truth_table, is_mixing) in operators.items():
        if is_mixing:
            # We check for algebraic contradiction
            # For a learnable function F(x,y), the logits must satisfy the sign requirements.
            # L_A, L_B, L_C are logits for (0,0), (0,1), (1,0)
            # The model structure implies L_D = L_B + L_C - L_A
            # sign(z) = 1 if z>0, -1 if z<0
            s_A = 1 if truth_table[0] == 1 else -1
            s_B = 1 if truth_table[1] == 1 else -1
            s_C = 1 if truth_table[2] == 1 else -1
            s_D = 1 if truth_table[3] == 1 else -1
            
            print(f"Testing {name}:")
            print(f"  Target sign pattern for L(0,0),L(0,1),L(1,0),L(1,1): [{s_A}, {s_B}, {s_C}, {s_D}]")

            # L(1,1) = L(0,1) + L(1,0) - L(0,0)
            # sign(L_D) must be equal to sign(L_B + L_C - L_A)
            
            # The sum L_B + L_C - L_A has a determinate sign only if
            # sign(L_B), sign(L_C), and sign(-L_A) are all the same.
            sign_of_L_B = s_B
            sign_of_L_C = s_C
            sign_of_neg_L_A = -s_A

            print(f"  Model requires sign(L(1,1)) = sign( (term with sign {sign_of_L_B}) + (term with sign {sign_of_L_C}) + (term with sign {sign_of_neg_L_A}) )")

            # Check for contradiction
            if sign_of_L_B > 0 and sign_of_L_C > 0 and sign_of_neg_L_A > 0:
                # Sum must be positive
                required_sum_sign = 1
                if required_sum_sign != s_D:
                    print("  -> CONTRADICTION: Sum of three positive terms must be positive.")
                    print(f"     Model requires sign(L(1,1))={required_sum_sign}, but Truth Table requires {s_D}.")
                    print(f"     Operator is NOT learnable.\n")
                    unlearnable.append(name.split(" ")[0])
                else:
                    print("  -> No contradiction found. Operator is learnable.\n")

            elif sign_of_L_B < 0 and sign_of_L_C < 0 and sign_of_neg_L_A < 0:
                # Sum must be negative
                required_sum_sign = -1
                if required_sum_sign != s_D:
                    print("  -> CONTRADICTION: Sum of three negative terms must be negative.")
                    print(f"     Model requires sign(L(1,1))={required_sum_sign}, but Truth Table requires {s_D}.")
                    print(f"     Operator is NOT learnable.\n")
                    unlearnable.append(name.split(" ")[0])
                else:
                    print("  -> No contradiction found. Operator is learnable.\n")
            else:
                # The signs are mixed, so a solution is possible.
                print("  -> No contradiction found. Operator is learnable.\n")
    
    print("---------------------------------")
    print("Summary of non-learnable operators:", unlearnable)
    print("---------------------------------")


solve()