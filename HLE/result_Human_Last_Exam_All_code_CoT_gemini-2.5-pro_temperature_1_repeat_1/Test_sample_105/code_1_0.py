import numpy as np

def run_analysis():
    """
    Analyzes which mixed-dimension logical operators can be learned by the given heuristic.
    """
    # Define the logical operators mixing dimensions of h1 and h2.
    # The operator h(p1, q2) takes the first dimension of h1 and the second of h2.
    # In our notation, p1=h1[0], q2=h2[1].
    operators = {
        "X' (p1 XOR q2)": lambda p1, q2: p1 ^ q2,
        "C' (p1 AND q2)": lambda p1, q2: p1 & q2,
        "D' (p1 OR q2)": lambda p1, q2: p1 | q2,
        "E' (p1 EQUIV q2)": lambda p1, q2: 1 if p1 == q2 else 0,
        "I' (p1 -> q2)": lambda p1, q2: 1 if not p1 or q2 else 0,
    }

    # The model's score function has the form: score = f(p1, p2) + g(q1, q2).
    # The target function is h(p1, q2). For the model to learn, the true/false
    # classes of h(p1, q2) must be linearly separable for ALL inputs (p1,q1,p2,q2).
    #
    # We test for contradiction. A contradiction occurs if:
    # sign(h(1,1) - h(1,0)) is opposite to sign(h(0,1) - h(0,0)).
    # This is because the first difference implies a constraint on g (e.g., g_for_q2=1 > g_for_q2=0)
    # and the second difference implies the opposite constraint on g (e.g., g_for_q2=1 < g_for_q2=0),
    # which is impossible to satisfy simultaneously with a single model.

    unlearnable_operators = []

    print("Analyzing which mixed-dimension operators cannot be learned.")
    print("-" * 65)
    print("A relation h(p1, q2) cannot be learned if the required constraints on")
    print("the learned weights lead to a logical contradiction.\n")

    for name, op in operators.items():
        # Evaluate the operator h for the four corners of the (p1, q2) input space.
        h_11 = op(1, 1)
        h_10 = op(1, 0)
        h_01 = op(0, 1)
        h_00 = op(0, 0)

        # The sign of this difference determines the required inequality for g when p1=1
        diff_when_p1_is_1 = h_11 - h_10
        # The sign of this difference determines the required inequality for g when p1=0
        diff_when_p1_is_0 = h_01 - h_00

        print(f"Operator: {name}")
        # Final Equation Part 1
        print(f"Analysis for p1=1: h(1,1) - h(1,0) = {h_11} - {h_10} = {diff_when_p1_is_1}")
        # Final Equation Part 2
        print(f"Analysis for p1=0: h(0,1) - h(0,0) = {h_01} - {h_00} = {diff_when_p1_is_0}")

        # Check for contradiction: signs are non-zero and opposite
        is_unlearnable = (diff_when_p1_is_1 * diff_when_p1_is_0) < 0

        if is_unlearnable:
            unlearnable_operators.append(name.split(" ")[0])
            print("Result: Contradiction found. The signs of the differences are opposite.")
            print(f"This operator CANNOT be learned.")
        else:
            print("Result: No contradiction found. This operator is learnable.")
        print("-" * 65)

    print("\nSummary:")
    print("All element-wise operators (X, C, D, E, I) are learnable.")
    print("The mixed-dimension operators that cannot be learned are:", unlearnable_operators)

run_analysis()