def solve_auditor_choice():
    """
    This function determines the set of firm types (theta) that choose lenient auditors
    under the scenario described in the problem.
    """

    # The problem describes a scenario where "companies keeping more truthful accounts
    # choose more lenient auditors".
    # In the context of theta={0: malpractice, 1: truthful} and x={0: lenient, 1: strict},
    # this implies a "reverse sorting" pattern:
    # - The truthful firm (theta=1) chooses the lenient auditor (x=0).
    # - The malpractice firm (theta=0) chooses the strict auditor (x=1).

    # Our derivation shows this scenario is only possible under one key condition derived from
    # the cost-minimization behavior of the malpractice firm (theta=0):
    # Condition for reverse sorting: F(0) < t(0)
    # where F(0) is the penalty for malpractice and t(0) is the tax liability for a malpractice firm.

    # Now, assuming this condition holds, we identify which firm types choose lenient auditors.

    set_of_lenient_choosers = set()

    # Analyze firm type theta = 0 (malpractice)
    # Under the condition F(0) < t(0), this firm's cost with a strict auditor is less than
    # with a lenient one. It will choose the strict auditor.
    # Therefore, theta=0 is not in the set.

    # Analyze firm type theta = 1 (truthful)
    # A truthful firm's cost is the same with either auditor type, C(1,0) = C(1,1).
    # This firm is indifferent. For the reverse sorting scenario to hold,
    # we must assume that an indifferent firm chooses the lenient auditor.
    # Therefore, theta=1 is in the set.
    set_of_lenient_choosers.add(1)

    print("The scenario implies a specific sorting pattern where more truthful firms choose more lenient auditors.")
    print("This requires the condition F(0) < t(0), meaning the penalty for malpractice is smaller than the tax liability.")
    print("Assuming this condition holds, we can determine the choices for each firm type:")
    print("- The firm with theta=0 (malpractice) chooses the strict auditor.")
    print("- The firm with theta=1 (truthful) is indifferent and is assumed to choose the lenient auditor.")
    print("\nTherefore, the set of values for theta where companies choose more lenient auditors is:")
    # Print each element of the set as requested
    for val in set_of_lenient_choosers:
        print(val)

solve_auditor_choice()