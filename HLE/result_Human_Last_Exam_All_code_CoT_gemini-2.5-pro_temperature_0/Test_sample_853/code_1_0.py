def analyze_corporate_structure():
    """
    Analyzes the given corporate structure options based on the partners' requirements.
    """
    requirements = {
        1: "Equal control over the company and its operations.",
        2: "Alex receives payments via salary.",
        3: "Tyler receives payments via dividends.",
        4: "Option for future non-voting investors who can receive a financial return."
    }

    analysis = {
        'A': "FAILS: The non-voting Class C shares are not eligible for dividends, making them unsuitable for future investors.",
        'B': "FAILS: Tyler is issued Class A shares which are not eligible for dividends, contradicting his requirement.",
        'C': "FAILS: The non-voting Class C shares offer no financial return (dividends/assets), making them unsuitable for future investors.",
        'D': "FAILS: This structure has no class of non-voting shares for future investors.",
        'E': "PASSES: This is the only option that satisfies all requirements. \n"
             "   - Control: Equal voting shares (50/50) provide equal fundamental control. \n"
             "   - Payments: Alex can take a salary, and Tyler's shares are dividend-eligible. \n"
             "   - Investors: It includes Class C non-voting, dividend-paying shares, which are perfect for future investors."
    }

    print("--- Analysis of Corporate Structure Options ---")
    for option, result in analysis.items():
        print(f"Option {option}: {result}\n")

    final_choice = 'E'
    print(f"Conclusion: Based on the analysis, Option {final_choice} is the only one that fulfills all four requirements.")

analyze_corporate_structure()
print("<<<E>>>")