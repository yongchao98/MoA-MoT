def solve_puzzle():
    """
    This function explains the logical steps to solve Fred's puzzle.
    """
    print("Analyzing the puzzle statements to find the single true statement:")
    print("-" * 50)

    # Statement 5 Analysis
    print("1. Analyzing S5: 'My best friend is older than me. He is 18 and I am 19.'")
    s5_is_false = True
    print("   - The statement asserts the friend is older, but gives ages (18 vs 19) that prove the friend is younger.")
    print("   - This is a self-contradiction, so the statement is logically FALSE.")
    print(f"   - Conclusion: The truth value of S5 is {not s5_is_false}.")
    print("-" * 50)

    # Statement 4 and 3 Analysis
    print("2. Analyzing S4 and S3:")
    print("   - S4: 'The number of the males is equal to number of the females in my class.'")
    print("   - S3: 'The total number of my friend in my class is even.'")
    print("   - If S4 is true (Males = Females), then the total number of friends is 2 * Males, which is always an even number.")
    print("   - This means that if S4 were true, S3 would also have to be true.")
    print("   - The puzzle rules state ONLY ONE statement is true. Therefore, S4 cannot be the true statement.")
    s4_is_false = True
    print(f"   - Conclusion: The truth value of S4 is {not s4_is_false}.")
    print("-" * 50)

    # Statement 1 Analysis
    print("3. Analyzing S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("   - This is a compound statement. Both parts must be true.")
    print("   - Part A: 'My name is Fred' -> True (from the problem setup).")
    print("   - Part B: 'if yesterday was after tomorrow, it would be Wednesday.'")
    print("   - The premise 'yesterday was after tomorrow' is impossible and therefore always false.")
    print("   - In logic, a conditional statement with a false premise is always true.")
    print("   - Since Part A and Part B are both true, the entire statement S1 is true.")
    s1_is_true = True
    print(f"   - Conclusion: The truth value of S1 is {s1_is_true}.")
    print("-" * 50)

    # Final Conclusion
    print("Summary of Truth Values:")
    print(f"S1: {s1_is_true}, S4: {not s4_is_false}, S5: {not s5_is_false}")
    print("Since only one statement can be true, it must be S1.")
    print("\nThe problem states Fred said the true statement on his truth-telling day.")
    print("The true statement, S1, is the only one that contains a day of the week: Wednesday.")
    print("\nTherefore, the day Fred tells the truth is Wednesday.")

solve_puzzle()
<<<C>>>