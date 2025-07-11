def solve_case():
    """
    Analyzes the corporate structure options to find the one that meets all four requirements.
    
    Requirements:
    1. Equal control for Alex and Tyler.
    2. Alex paid by salary.
    3. Tyler paid by dividends.
    4. Option for non-voting investors.

    Analysis of Options:
    - A: Fails on requirement 4. The non-voting Class C shares are worthless and not viable for investors.
    - B: Fails on requirement 3. Tyler is given Class A shares which are not eligible for dividends.
    - C: Fails on requirement 4. The non-voting Class C shares are worthless.
    - D: Fails on requirement 4. It does not authorize any class of non-voting shares for future investors.
    - E: The share structure itself is perfect. It provides voting shares (Class A/B) for the founders and non-voting, dividend-paying shares (Class C) for future investors. This allows all four requirements to be met. The initial appointment of extra directors is an operational issue that the founders could correct with their equal voting power, it is not a flaw in the fundamental share structure.
    """
    # The correct option is E because its share structure is the only one that facilitates all four requirements.
    correct_answer = 'E'
    
    print(f"The analysis of the requirements leads to the conclusion that Option E provides the most suitable corporate structure.")
    print(f"1) Equal Control: Possible with equal holdings of voting Class A and Class B shares.")
    print(f"2) Alex's Salary: Possible as a director/employee.")
    print(f"3) Tyler's Dividends: Possible as a holder of dividend-paying Class A or B shares.")
    print(f"4) Non-voting Investors: Possible through the issuance of non-voting, dividend-paying Class C shares.")
    print(f"Therefore, the correct answer is {correct_answer}.")

solve_case()
<<<E>>>