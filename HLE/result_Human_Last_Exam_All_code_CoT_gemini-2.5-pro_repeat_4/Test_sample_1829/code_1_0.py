def solve_puzzle():
    """
    This script solves the logic puzzle by analyzing each statement step-by-step.
    """

    print("Step 1: Analyzing the truth value of each statement.")
    print("--------------------------------------------------")

    # Statement S5 Analysis
    print("Analyzing S5: 'My best friend is older than me. He is 18 and I am 19.'")
    my_age = 19
    friend_age = 18
    is_friend_older = friend_age > my_age
    print(f"I am {my_age} and my friend is {friend_age}. The claim that my friend is older is {is_friend_older}.")
    print("S5 is self-contradictory, so it must be FALSE.\n")

    # Statement S4 and S3 Analysis
    print("Analyzing S4 and S3 together.")
    print("S4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("If S4 is true, then Total_Friends = Males + Females = 2 * Males, which is always an even number.")
    print("This means if S4 is true, S3 must also be true.")
    print("The puzzle states only ONE statement is true. If S4 were true, we would have at least two true statements (S3 and S4).")
    print("Therefore, S4 must be FALSE.\n")

    # Statement S1 Analysis
    print("Analyzing S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("This is a conditional statement. Let's look at the premise: 'yesterday was after tomorrow'.")
    print("Interpreting this alphabetically, the word 'yesterday' comes after 'tomorrow' in the dictionary.")
    premise_is_true = "yesterday" > "tomorrow"
    print(f"Is the premise 'yesterday' > 'tomorrow' true? {premise_is_true}.")
    print("So, the premise of the conditional is TRUE.")
    print("For a conditional statement ('if P, then Q') to be true when the premise (P) is true, the conclusion (Q) must also be true.")
    print("The conclusion is 'it would be Wednesday'.")
    print("Therefore, for S1 to be true, it must be spoken on a Wednesday. On any other day, S1 is false.\n")

    print("Step 2: Identifying the true statement and the day.")
    print("--------------------------------------------------")
    print("We have determined that S4 and S5 are false.")
    print("We are told exactly one statement is true, so that statement must be S1.")
    print("Since S1 is only true when spoken on a Wednesday, the day of the statement must be WEDNESDAY.\n")

    print("Step 3: Confirming with the final condition.")
    print("---------------------------------------------")
    print("The puzzle states: 'Only one statement is true and Fred said it on the day when he tells the truth.'")
    print("We found the true statement is S1, and it was said on Wednesday.")
    print("This means Wednesday must be a day when Fred tells the truth.\n")
    
    print("Final Answer Derivation:")
    print("The single true statement is S1.")
    print("S1 is only true on Wednesday.")
    print("Fred uttered the true statement on his truth-telling day.")
    print("Therefore, the day Fred tells the truth (in this context) is Wednesday.")
    
solve_puzzle()
<<<C>>>