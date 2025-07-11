def solve_puzzle():
    """
    This function solves the logic puzzle by analyzing each statement
    and printing the step-by-step reasoning.
    """

    print("Step 1: Analyzing the statements to determine their logical truth value.")
    print("-" * 60)

    # Analysis of S5
    print("Analyzing S5: 'My best friend is older than me. He is 18 and I am 19.'")
    print("A person who is 18 years old is younger than a person who is 19 years old.")
    print("The statement claims the 18-year-old is older, which is a direct contradiction.")
    print("Therefore, statement S5 is logically FALSE.")
    print("-" * 60)

    # Analysis of S4 and S3
    print("Analyzing S4: 'The number of the males is equal to number of the females in my class.'")
    print("And S3: 'The total number of my friend in my class is even.'")
    print("If S4 were true, the total number of friends would be (number of males + number of females).")
    print("Since males = females, the total would be (males + males) = 2 * males, which is always an even number.")
    print("This means that if S4 is true, S3 must also be true.")
    print("However, the puzzle states that ONLY ONE statement is true.")
    print("Therefore, it's impossible for S4 to be true.")
    print("So, statement S4 must be FALSE.")
    print("-" * 60)

    # Analysis of S1
    print("Analyzing S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("This statement is a compound statement: (Part A) AND (Part B).")
    print("Part A is 'My name is Fred', which we can assume is true.")
    print("Part B is 'if yesterday was after tomorrow, it would be Wednesday.'")
    print("The condition 'yesterday was after tomorrow' is a logical impossibility; it can never happen. Thus, it is always FALSE.")
    print("In formal logic, an 'if-then' statement with a false condition (e.g., 'if FALSE, then ...') is always TRUE, regardless of the outcome.")
    print("So, Part B is logically TRUE.")
    print("Since Part A (True) and Part B (True) are both true, the entire statement S1 is TRUE.")
    print("-" * 60)

    # Consolidating the results
    print("Step 2: Consolidating the truth values.")
    print("We have determined the following:")
    print("S1 is TRUE.")
    print("S5 is FALSE.")
    print("S4 is FALSE.")
    print("Since only one statement can be true, S3 must also be FALSE.")
    print("\nSummary of Truth Values:")
    print("S1: TRUE")
    print("S3: FALSE")
    print("S4: FALSE")
    print("S5: FALSE")
    print("This aligns with the puzzle's rule that exactly one statement is true.")
    print("-" * 60)

    # Determining the day
    print("Step 3: Determining Fred's truth-telling day.")
    print("The puzzle states that Fred said the one true statement on the day he tells the truth.")
    print("The single true statement is S1.")
    print("S1 is: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("In this type of logic puzzle, the specific information contained within the key statement provides the answer.")
    print("The only information in the true statement that refers to a day of the week is 'Wednesday'.")
    print("\nConclusion: The day Fred tells the truth is Wednesday.")

solve_puzzle()
<<<C>>>