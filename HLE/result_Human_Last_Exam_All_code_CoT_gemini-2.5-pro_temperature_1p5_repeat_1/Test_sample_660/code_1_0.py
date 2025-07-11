def solve():
    """
    Solves for the set of theta values based on the problem description.
    """
    # Define symbolic variables for the parameters
    t0 = "t(0)"
    p1 = "p(1)"
    F0 = "F(0)"

    print("To solve the problem, we establish the conditions under which a 'sorting equilibrium' exists, where truthful firms choose lenient auditors and malpractice firms choose strict auditors. This matches the problem's premise that 'companies keeping more truthful accounts choose more lenient auditors'.")
    print("-" * 80)

    # Malpractice Firm's choice
    print("1. Analyze the choice for a malpractice firm (theta = 0):")
    print(f"   - The cost of hiring a lenient auditor (x=0) is certain, as this auditor always issues a report. The cost is the tax liability:")
    print(f"     Cost(theta=0, x=0) = {t0}\n")

    print(f"   - The cost of hiring a strict auditor (x=1) is the expected cost after the auditor refuses the report (which they always do for theta=0). This cost depends on the probability of a Revenue Service audit, p(1).")
    print(f"     Cost(theta=0, x=1) = {p1} * ({t0} + {F0})\n")

    print("   - For a malpractice firm to choose a strict auditor, the cost must be lower. This gives us the condition for the sorting equilibrium:")
    print(f"     Condition: {p1} * ({t0} + {F0}) < {t0}")
    print("-" * 80)

    # Truthful Firm's choice
    print("2. Analyze the choice for a truthful firm (theta = 1):")
    print("   - A truthful firm will always have its report issued by any auditor, as the probability of refusal is x*(1-1) = 0.")
    print("   - Therefore, its cost is always t(1), and it is indifferent between auditors. We assume it chooses the lenient auditor.")
    print("-" * 80)
    
    # Conclusion
    print("3. Conclusion:")
    print("   Assuming the condition for the sorting equilibrium holds, we can identify the set of firms that choose lenient auditors:")
    print("   - Firms with theta=1 choose lenient auditors.")
    print("   - Firms with theta=0 choose strict auditors.")
    print("\nTherefore, the set of values of theta for which companies choose lenient auditors is:")
    
    final_set = "{1}"
    print(final_set)
    

solve()
