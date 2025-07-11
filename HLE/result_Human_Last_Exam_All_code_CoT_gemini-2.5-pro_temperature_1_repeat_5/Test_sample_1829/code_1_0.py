def solve_puzzle():
    """
    This function solves the logic puzzle by analyzing each statement.
    """
    print("Let's analyze the statements to find the single true one.")

    # --- Statement 5 Analysis ---
    # This is the most straightforward statement to evaluate.
    my_age = 19
    friend_age = 18
    print("\n--- Analysis of S5: 'My best friend is older than me. He is 18 and I am 19.' ---")
    print(f"This statement claims the friend is older, but gives the ages as Fred = {my_age} and Friend = {friend_age}.")
    # The puzzle requires outputting the numbers in the equation.
    print(f"The factual comparison is: {my_age} > {friend_age}.")
    print("This contradicts the claim that the friend is older. So, S5 is FALSE.")

    # --- Statement 1 Analysis ---
    # This statement requires understanding a bit of formal logic.
    print("\n--- Analysis of S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.' ---")
    print("This is a compound statement. Let's look at the second part: 'if yesterday was after tomorrow...'")
    print("The condition 'yesterday was after tomorrow' is a logical impossibility; it can never be true. It is FALSE.")
    print("In logic, an 'if-then' statement with a false condition is always considered TRUE, regardless of the outcome.")
    print("Since the first part ('My name is Fred') is true and the second part is logically true, the entire statement S1 is TRUE.")

    # --- Statements 3 & 4 Analysis ---
    # These statements are linked.
    print("\n--- Analysis of S3 & S4 ---")
    print("S4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("If S4 were true (males == females), the total number of friends would be (males + males) = 2 * males, which is always an EVEN number.")
    print("This would force S3 to also be true.")
    print("However, the puzzle states only ONE statement is true. Therefore, S4 cannot be true.")
    print("Since we know S1 is true, both S3 and S4 must be FALSE.")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print("Our analysis shows:")
    print("S1 is TRUE.")
    print("S3 is FALSE.")
    print("S4 is FALSE.")
    print("S5 is FALSE.")
    print("\nThe single true statement is S1.")
    print("The puzzle states Fred tells the truth on the day he makes his true statement.")
    print("The true statement, S1, contains the name of a day: 'Wednesday'.")
    print("\nTherefore, the day Fred tells the truth is Wednesday.")

solve_puzzle()
<<<C>>>