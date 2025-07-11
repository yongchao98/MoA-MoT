def solve_puzzle():
    """
    This function prints the step-by-step logical deduction to solve the puzzle.
    """
    print("Analyzing the statements to find the single true one:\n")

    # --- Step 1: Analyze Statement S5 ---
    print("Step 1: Analyzing S5")
    print("Statement S5: 'My best friend is older than me. He is 18 and I am 19.'")
    friend_age = 18
    my_age = 19
    print(f"The statement provides the ages: Best friend = {friend_age}, Me (Fred) = {my_age}.")
    print(f"The comparison shows that {friend_age} < {my_age}, which means the best friend is younger.")
    print("This contradicts the first part of the statement ('My best friend is older than me').")
    print("Conclusion: Statement S5 is internally contradictory and therefore must be FALSE.\n")

    # --- Step 2: Analyze Statements S3 and S4 ---
    print("Step 2: Analyzing S3 and S4")
    print("Statement S4: 'The number of the males is equal to number of the females in my class.'")
    print("Statement S3: 'The total number of my friend in my class is even.'")
    print("Let M be the number of males and F be the number of females.")
    print("If S4 is true, then M = F.")
    print("If M = F, the total number of friends is M + F = M + M = 2 * M.")
    print("A number of the form 2 * M is always even.")
    print("This means that if S4 is TRUE, then S3 must also be TRUE.")
    print("However, the puzzle states that ONLY ONE statement is true.")
    print("Therefore, the condition 'S4 is true' is impossible as it leads to two true statements.")
    print("Conclusion: Statement S4 must be FALSE.\n")

    # --- Step 3: Analyze Statement S1 ---
    print("Step 3: Analyzing S1")
    print("Statement S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("This is a compound statement. Let's analyze its two parts.")
    print("Part A: 'My name is Fred.' The puzzle starts with 'Fred tells lies...', so we accept his name is Fred. This part is TRUE.")
    print("Part B: 'if yesterday was after tomorrow, it would be Wednesday.' This is a logical implication (if P then Q).")
    print("   - The premise P is 'yesterday was after tomorrow'. This is impossible and thus logically FALSE.")
    print("   - In logic, an 'if-then' statement with a false premise is always considered TRUE, regardless of the conclusion.")
    print("   - Therefore, Part B is TRUE.")
    print("Since S1 is (Part A AND Part B), and both parts are TRUE, the entire statement is TRUE.")
    print("Conclusion: Statement S1 is logically TRUE.\n")

    # --- Step 4: Final Conclusion ---
    print("Step 4: Synthesizing and Final Conclusion")
    print("We have determined the following:")
    print(" - S1 is TRUE.")
    print(" - S4 is FALSE.")
    print(" - S5 is FALSE.")
    print("The puzzle requires exactly one true statement. Our analysis shows S1 is that statement.")
    print("(This also means S3 must be false).")
    print("\nThe single true statement, which Fred spoke on his truth-telling day, is S1.")
    print("The statement S1 mentions the day 'Wednesday'.")
    print("In this classic puzzle format, the content of the true statement reveals the answer.")
    print("\nTherefore, the day Fred tells the truth is Wednesday.")

solve_puzzle()
<<<C>>>