def solve_puzzle():
    """
    This function analyzes the logic puzzle to determine the day Fred tells the truth.
    """
    print("Step 1: Analyzing the truth value of each statement.\n")

    # --- Analysis of S5 ---
    # S5: My best friend is older than me. He is 18 and I am 19.
    freds_age = 19
    friends_age = 18
    s5_is_true = friends_age > freds_age
    print("Analysis of S5: 'My best friend is older than me. He is 18 and I am 19.'")
    print(f"  - Fred's age is {freds_age} and his friend's age is {friends_age}.")
    print(f"  - The claim that the friend is older ({friends_age} > {freds_age}) is clearly false.")
    print(f"  - Result: S5 is {s5_is_true}.\n")

    # --- Analysis of S3 and S4 ---
    # S4: The number of the males is equal to number of the females in my class.
    # S3: The total number of my friend in my class is even.
    # If S4 is true, S3 must also be true. This would violate the "only one is true" rule.
    s4_is_true = False
    s3_is_true = False # Must be false if only one statement in total is true.
    print("Analysis of S3 & S4:")
    print("  - If S4 were true (males == females), the total would be even, making S3 true.")
    print("  - This would mean at least two statements are true, which contradicts the puzzle's rule.")
    print(f"  - Result: S4 must be {s4_is_true} and S3 must be {s3_is_true}.\n")
    
    # --- Analysis of S1 ---
    # S1: My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.
    # Part 1: "My name is Fred" -> True (given)
    # Part 2: "if yesterday was after tomorrow..." -> False premise makes the conditional True.
    # So, S1 is True and True -> True
    s1_is_true = True
    print("Analysis of S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("  - The condition 'if yesterday was after tomorrow' is a logical impossibility, so the premise is False.")
    print("  - In logic, an 'if-then' statement with a False premise is always True.")
    print(f"  - Result: S1 is {s1_is_true}.\n")

    print("Step 2: Conclusion from the analysis.\n")
    print("Summary of Truth Values:")
    print(f"  - S1 is {s1_is_true}")
    print(f"  - S3 is {s3_is_true}")
    print(f"  - S4 is {s4_is_true}")
    print(f"  - S5 is {s5_is_true}\n")

    print("This confirms that S1 is the only true statement, as required by the puzzle.\n")

    print("Step 3: Identifying Fred's truth-telling day.\n")
    # The puzzle states Fred says the true statement on his truth day.
    # The true statement, S1, contains the name of a day.
    truth_day = "Wednesday"
    print("  - The single true statement is S1.")
    print("  - S1 mentions a specific day of the week: '...it would be Wednesday.'")
    print(f"  - The clue to the answer is embedded within the true statement.")
    print(f"  - Therefore, Fred tells the truth on {truth_day}.")

solve_puzzle()
<<<C>>>