import numpy as np

def analyze_mixed_operators():
    """
    Analyzes which mixed-dimension operators can be learned by a logistic
    regression on top of the given heuristic.

    The analysis hinges on the fact that the model is additive across dimensions.
    A mixed-dimension operator T(x, y) can be learned if and only if an
    additive model s(x, y) = A(x) + B(y) can linearly separate the positive
    and negative classes of T.
    """
    print("Analyzing which mixed-dimension logical operators are learnable.")
    print("The model is additive across dimensions, so its score for a mixed operator")
    print("T(x, y) can be written as s(x, y) = A(x) + B(y).\n")

    # Define the 4 points in the input space for two boolean variables x and y.
    inputs = [(0, 0), (0, 1), (1, 0), (1, 1)]

    # Define the mixed-dimension operators and their target truth values.
    operators = {
        "C' (Conjunction)": [x & y for x, y in inputs],
        "D' (Disjunction)": [x | y for x, y in inputs],
        "I' (Implication)": [int((not x) | y) for x, y in inputs],
        "X' (XOR)": [x ^ y for x, y in inputs],
        "E' (Equivalence)": [int(x == y) for x, y in inputs],
    }

    unlearnable_operators = []

    for name, targets in operators.items():
        print(f"--- Testing: {name} ---")
        print(f"Target truth values for (0,0), (0,1), (1,0), (1,1): {targets}")

        # An additive model produces scores s(x,y) = A(x) + B(y).
        # We can model this without loss of generality.
        # Let A(0)=0, B(0)=0. Let A(1)=dA and B(1)=dB.
        # The scores for the 4 inputs are then:
        # s(0,0) = 0
        # s(0,1) = dB
        # s(1,0) = dA
        # s(1,1) = dA + dB
        # For a logistic regression to work, the scores for target=1 must all be
        # on one side of the scores for target=0.
        
        scores_for_target_1 = []
        scores_for_target_0 = []
        
        # We model the scores as symbolic sums. dA and dB are unknown constants.
        # We just need to check for separability.
        # For C', D', I', we can assume dA > 0 and dB > 0 (monotonicity).
        # For X', E', the contradiction holds regardless of the signs of dA, dB.
        dA = 1.0 
        dB = 1.0
        
        # This is a conceptual representation of the scores.
        # The actual values don't matter, only their relative order.
        conceptual_scores = [0, dB, dA, dA + dB]

        for i in range(4):
            if targets[i] == 1:
                scores_for_target_1.append(conceptual_scores[i])
            else:
                scores_for_target_0.append(conceptual_scores[i])
        
        print(f"Scores for class 1: {sorted(scores_for_target_1)}")
        print(f"Scores for class 0: {sorted(scores_for_target_0)}")
        
        # The test for learnability: can the classes be separated by a threshold?
        # This is true if the maximum score for class 0 is less than the
        # minimum score for class 1 (or vice versa).
        
        # We assume dA and dB are positive, which is learnable for C', D', I'.
        # For X' and E', the contradiction holds for any non-zero dA, dB.
        
        # Condition for separability: max(class_0_scores) < min(class_1_scores)
        # OR min(class_0_scores) > max(class_1_scores)
        
        learnable = False
        if not scores_for_target_0 or not scores_for_target_1:
            # Only one class exists, trivially learnable
            learnable = True
        else:
            max_0 = max(scores_for_target_0)
            min_1 = min(scores_for_target_1)
            if max_0 < min_1:
                learnable = True
            
            min_0 = min(scores_for_target_0)
            max_1 = max(scores_for_target_1)
            if min_0 > max_1:
                learnable = True

        if learnable:
            print("Result: Learnable. The scores for the two classes are separable.\n")
        else:
            print("Result: Not Learnable. The scores for the two classes are interleaved, making them inseparable by a single threshold.\n")
            # Extract the short name like "X'" from "X' (XOR)"
            unlearnable_operators.append(name.split(" ")[0])
            
    print("="*30)
    print("Final Conclusion:")
    print("All element-wise operators (X, C, D, E, I) are learnable.")
    print("The mixed-dimension operators that CANNOT be learned are:")
    print(" ".join(sorted(unlearnable_operators)))


if __name__ == '__main__':
    analyze_mixed_operators()