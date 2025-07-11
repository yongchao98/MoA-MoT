def solve_fred_puzzle():
    """
    Solves the logic puzzle by analyzing each statement to find the single true one.
    """
    print("Analyzing the statements to determine which one is true:\n")

    # --- S5 Analysis ---
    # Statement: "My best friend is older than me. He is 18 and I am 19."
    my_age = 19
    friend_age = 18
    # The statement claims the friend is older, but the ages contradict this.
    # We can model the claim as an inequality: friend_age > my_age
    is_friend_older = friend_age > my_age
    print("Analysis of S5: 'My best friend is older than me. He is 18 and I am 19.'")
    print(f"Evaluating the claim that the friend is older using the numbers provided:")
    print(f"Is {friend_age} > {my_age}? This evaluates to {is_friend_older}.")
    print("Result: Statement S5 contradicts itself, so it is FALSE.\n")

    # --- S4 and S3 Analysis ---
    print("Analysis of S3 & S4:")
    print("S4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("Logic: If S4 is true, then Total Friends = (Number of Males) + (Number of Females). If we assume Number of Males = X, then Total Friends = X + X = 2 * X.")
    print("A number multiplied by 2 is always even. Therefore, if S4 is true, S3 must also be true.")
    print("However, the puzzle states that only ONE statement is true. Having both S3 and S4 be true is a contradiction.")
    print("Result: S4 cannot be the true statement, so it must be FALSE.\n")

    # --- S1 Analysis ---
    print("Analysis of S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("The condition 'if yesterday was after tomorrow...' relies on a premise that is impossible. It is a logical falsehood.")
    print("In logic, an 'if...then' statement with a premise that is always false is known as a 'vacuously true' statement.")
    print("Result: Statement S1 is logically TRUE, regardless of what follows the 'then'.\n")

    # --- Final Conclusion ---
    print("Conclusion:")
    print("S1 is TRUE, while S4 and S5 are provably FALSE.")
    print("Since only one statement can be true, it must be S1.")
    print("The true statement, S1, mentions the day 'Wednesday'. This is the only clue to the specific day.")
    print("\nTherefore, the day Fred tells the truth is Wednesday.")

solve_fred_puzzle()
<<<C>>>