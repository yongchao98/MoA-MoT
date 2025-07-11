def solve_freds_puzzle():
    """
    This program analyzes the logic puzzle about Fred to determine the day
    he tells the truth. It explains the reasoning step by step.
    """
    print("Step 1: Analyzing the statements based on the puzzle's rules.")
    print("Rule: Exactly one of the four statements is true.\n")

    # --- Analysis of Statement 5 ---
    print("Analysis of S5: 'My best friend is older than me. He is 18 and I am 19.'")
    my_age = 19
    friend_age = 18
    # The core claim is that the friend is older, which we can write as an equation.
    is_friend_older = friend_age > my_age
    print(f"The statement includes two numbers: Fred's age ({my_age}) and his friend's age ({friend_age}).")
    print(f"The claim 'my best friend is older than me' can be written as the equation: {friend_age} > {my_age}")
    print(f"This equation evaluates to {is_friend_older}.")
    print("Because the statement claims something that is contradicted by the numbers it provides, S5 is logically FALSE.\n")

    # --- Analysis of Statements 3 and 4 ---
    print("Analysis of S3 and S4:")
    print("S4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("If S4 were true, the number of males would equal females (M=F).")
    print("In that case, the total number of friends would be M + F = 2*M, which is always an even number.")
    print("So, if S4 is true, S3 must also be true.")
    print("This would mean two statements are true, violating the rule that 'only one statement is true'.")
    print("Therefore, S4 must be FALSE.\n")

    # --- Analysis of Statement 1 ---
    print("Analysis of S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("The condition 'yesterday was after tomorrow' is impossible.")
    print("In puzzles like this, a statement 'if [impossible], then [X]' is a cryptic way of saying 'X is true'.")
    print("So, we interpret S1 as a claim that 'Today is Wednesday'.")
    print("This means S1 is TRUE if and only if today is Wednesday.\n")

    # --- Final Deduction ---
    print("Step 2: Identifying the single true statement and the day.")
    print("We have proven that S4 and S5 are FALSE.")
    print("The one true statement must be either S1 or S3.")
    print("If S3 were the true statement, we would have no way to identify the specific day.")
    print("For the puzzle to have a unique solution, S1 must be the true statement.")
    print("If S1 is true, then our interpretation means the day must be Wednesday.")
    print("\nConclusion: The single true statement is S1, which Fred said on Wednesday.")

    final_answer_day = "Wednesday"
    final_answer_choice = "C"

    print(f"\nThe day Fred said the true statement is {final_answer_day}.")

solve_freds_puzzle()
<<<C>>>