def solve_puzzle():
    """
    This function solves the logic puzzle by analyzing the statements step-by-step.
    """
    print("Step 1: Analyzing the statements to find the single true statement.")
    
    # Analysis of S5
    s5_is_false = True
    print("S5: 'My best friend is older than me. He is 18 and I am 19.'")
    print("This is a self-contradiction. A person who is 18 is younger than a person who is 19.")
    print("Therefore, S5 is FALSE.")
    print("-" * 20)

    # Analysis of S3 and S4
    print("S4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("If S4 is true, let the number of males be M and females be F. So, M = F.")
    print("The total number of friends would be M + F = M + M = 2 * M.")
    print("2 * M is always an even number. So, if S4 is true, S3 must also be true.")
    s4_implies_s3 = True
    
    print("The puzzle states ONLY ONE statement is true.")
    print("If S4 were the true statement, S3 would also be true, which is a contradiction (2 true statements).")
    s4_is_false = True
    print("Therefore, S4 must be FALSE.")
    print("-" * 20)

    # Conclusion on which statement is true
    print("S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("The condition 'yesterday was after tomorrow' is a logical impossibility, so it's always false.")
    print("In logic, an 'if P then Q' statement is always true if the premise P is false.")
    print("So, the second part of S1 is always true. The first part, 'My name is Fred', is a given.")
    print("This makes the entire statement S1 logically TRUE.")
    s1_is_true = True
    
    print("\nSince S4 and S5 are false, and only one statement can be true, S1 must be the one true statement.")
    print("This means S3 must also be false.")
    print("\nSummary of Truth Values:")
    print("S1: TRUE")
    print("S3: FALSE")
    print("S4: FALSE")
    print("S5: FALSE")
    print("-" * 20)
    
    print("Step 2: Interpreting S1 to find the day.")
    print("The question asks for a single day of the week, which implies Fred has only ONE truth-telling day.")
    print("The puzzle must be solved by finding the hidden meaning in the true statement, S1.")
    print("In logic puzzles, an impossible premise like 'if yesterday was after tomorrow...' often refers to the speaker's special or alternate state.")
    print("Since Fred lies on 6 days and tells the truth on 1 day, his 'special day' is his truth day.")
    print("So, we can decode S1 as: 'My truth day is Wednesday.'")
    print("-" * 20)

    print("Step 3: Determining the truth day.")
    print("We know that the decoded S1, 'My truth day is Wednesday', is the one true statement.")
    print("For this statement to be true, Fred's truth day must be Wednesday.")
    print("The puzzle also says he uttered the true statement on his truth day. This is consistent: on Wednesday, he truthfully states that his truth day is Wednesday.")
    final_day = "Wednesday"
    print(f"\nConclusion: The day Fred tells the truth is {final_day}.")

solve_puzzle()
<<<C>>>