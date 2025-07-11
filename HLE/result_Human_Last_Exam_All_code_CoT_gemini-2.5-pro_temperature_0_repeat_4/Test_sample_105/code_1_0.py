import numpy as np

def check_learnability(op_name, targets):
    """
    Checks if an operator is learnable by the additive model.

    The model's score L(x_i, y_j) has the structure L(1,1) = L(1,0) + L(0,1) - L(0,0).
    Let L(0,0) = 0, L(1,0) = A, L(0,1) = B. Then L(1,1) = A + B.
    The four points in the logit space (relative to L(0,0)) are [0, B, A, A+B].
    The problem is learnable if these four points are linearly separable according to the target labels.
    This is the classic test for linear separability of a 2-input boolean function.
    A 2-input boolean function is linearly separable iff:
    (t0 AND t1 AND t2 AND t3) is false, AND
    (t0 OR t1 OR t2 OR t3) is true, AND
    (t0 == t3) == (t1 == t2)
    where t0,t1,t2,t3 are the target values for (0,0), (0,1), (1,0), (1,1).
    """
    t0, t1, t2, t3 = targets[0], targets[1], targets[2], targets[3]

    # Condition 1: Not all targets are the same
    all_same = (t0 == t1 and t1 == t2 and t2 == t3)
    if all_same:
        # Trivially learnable, but not a useful function.
        # Let's assume the relations are non-trivial.
        is_learnable = True
    else:
        # Condition for linear separability of a 2-input boolean function
        is_learnable = (t0 == t3) == (t1 == t2)

    print(f"Operator: {op_name}")
    print(f"  Targets for (x_i,y_j) in [(0,0), (0,1), (1,0), (1,1)]: {targets}")
    if is_learnable:
        print("  Result: Learnable")
    else:
        print("  Result: NOT Learnable")
    print("-" * 20)
    return is_learnable

def main():
    # Define the target outputs for each mixing-dimension operator
    # Targets are for (x_i, y_j) inputs: (0,0), (0,1), (1,0), (1,1)
    operators = {
        "X' (XOR mixing)":          [0, 1, 1, 0],
        "C' (Conjunction mixing)":  [0, 0, 0, 1],
        "D' (Disjunction mixing)":  [0, 1, 1, 1],
        "E' (Equivalence mixing)":  [1, 0, 0, 1],
        "I' (Implication mixing)":  [1, 1, 0, 1],
    }

    not_learnable_ops = []
    for name, targets in operators.items():
        if not check_learnability(name, targets):
            # Extract the short name like X'
            short_name = name.split(" ")[0]
            not_learnable_ops.append(short_name)

    print("\nSummary:")
    print("The operators that CANNOT be learned are: " + ", ".join(not_learnable_ops))

if __name__ == "__main__":
    main()