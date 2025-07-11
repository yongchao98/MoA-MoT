def analyze_corporate_structure():
    """
    Analyzes the provided options to find the best corporate structure.
    """
    print("Analyzing the requirements for the corporate structure:")
    print("1. Equal Control for Alex and Tyler.")
    print("2. Alex receives a salary.")
    print("3. Tyler receives dividends.")
    print("4. Future investors can be brought in without voting rights.\n")

    print("Evaluating the options:")
    print("Option A: Fails. The non-voting Class C shares have no value (no dividends, no assets) and cannot attract investors.")
    print("Option B: Fails. Tyler is issued Class A shares, which are not eligible for dividends, contradicting a key requirement.")
    print("Option C: Fails. The non-voting Class C shares have no value, making them unsuitable for future investors.")
    print("Option D: Fails. It lacks a class of non-voting shares, so there is no way to bring in investors without giving them control.")
    print("Option E: Succeeds. It is the only structure that meets all needs:")
    print("  - Alex and Tyler receive voting shares (Class A & B). They can pay Tyler dividends and Alex a salary.")
    print("  - It includes a Class C share class that is non-voting but eligible for dividends, which is perfect for future investors.")
    print("  - Although the initial board is unequal, Alex and Tyler have equal voting shares (50 each), which gives them equal ultimate control to change the board.")

    final_answer = 'E'
    print(f"\nThe best corporate structure that satisfies all requirements is Option E.")

analyze_corporate_structure()
<<<E>>>