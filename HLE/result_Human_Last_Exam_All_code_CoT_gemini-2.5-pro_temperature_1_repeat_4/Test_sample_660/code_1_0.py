def solve_auditor_choice():
    """
    Solves for the set of theta for which more truthful firms
    choose more lenient auditors.
    """
    # Define the parameters of the model.
    # These values are chosen to be consistent with the problem statement's
    # requirement that t, p, and F are decreasing functions.
    # t(theta): Tax liability
    t = {0: 100, 1: 50}
    # p(x): Probability of Revenue Service audit
    p = {0: 0.5, 1: 0.2}
    # F(theta): Penalty for malpractice
    F = {0: 300, 1: 10}

    print("This script determines the set of firm types (theta) for which truthful firms (theta=1) choose lenient auditors (x=0) and malpracticing firms (theta=0) choose strict auditors (x=1).\n")
    print("Model Parameters:")
    print(f"  t(0) = {t[0]}, t(1) = {t[1]}")
    print(f"  p(1) = {p[1]}")
    print(f"  F(0) = {F[0]}\n")

    # The condition for a malpracticing firm (theta=0) to choose a strict auditor (x=1)
    # is C(0,1) < C(0,0), which simplifies to p(1)*F(0) < t(0)*(1-p(1)).
    # A truthful firm (theta=1) is indifferent and is assumed to choose the lenient auditor.
    
    # Let's check the inequality:
    lhs = p[1] * F[0]
    rhs = t[0] * (1 - p[1])

    print("Checking the condition for malpracticing firms to choose strict auditors:")
    print(f"  p(1) * F(0)  <  t(0) * (1 - p(1))")
    print(f"  {p[1]} * {F[0]}  <  {t[0]} * (1 - {p[1]})")
    print(f"  {lhs}  <  {t[0]} * {1 - p[1]}")
    print(f"  {lhs}  <  {rhs}\n")

    if lhs < rhs:
        result_set = "{0, 1}"
        print(f"The condition is TRUE.")
        print("This means malpracticing firms choose strict auditors, while truthful firms choose lenient ones.")
        print(f"The set of theta values for which this sorting occurs is {result_set}.")
        final_answer = result_set
    else:
        result_set = "{}" # Empty set
        print(f"The condition is FALSE.")
        print("This means malpracticing firms also choose lenient auditors (or are indifferent).")
        print(f"The described sorting pattern does not occur. The resulting set is the empty set: {result_set}.")
        final_answer = result_set
    
    # The final answer is the set determined by the inequality.
    # For this specific numerical example, we determine the outcome.
    # The final line <<<...>>> should contain the answer for the given parameters.
    print(f"\n<<<Result for the given parameters>>>\n{final_answer}")


solve_auditor_choice()