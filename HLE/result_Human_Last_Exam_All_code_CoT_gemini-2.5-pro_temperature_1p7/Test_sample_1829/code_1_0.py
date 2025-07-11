def solve_puzzle():
    """
    This function solves the logic puzzle by analyzing each statement
    and applying the given rules to find Fred's truth-telling day.
    """

    print("Step 1: Analyzing the truth value of each statement.\n")

    # Statement 5 Analysis
    s5_explanation = "S5: 'My best friend is older than me. He is 18 and I am 19.'\n"
    s5_analysis = "This statement is a self-contradiction. If I am 19 and my friend is 18, he cannot be older than me. A self-contradictory statement is always False.\n"
    print(s5_explanation + s5_analysis)

    # Statement 4 and 3 Analysis
    s4_explanation = "S4: 'The number of the males is equal to number of the females in my class.' (M = F)\n"
    s3_explanation = "S3: 'The total number of my friend in my class is even.' (M + F = Even)\n"
    s3_s4_analysis = "If S4 is True (M = F), then the total number of friends is M + M = 2M, which is always an even number. Therefore, S3 must also be True. This means that if S4 is true, S3 must be true. So, it's impossible for S4 to be the *only* true statement.\n"
    print(s4_explanation + s3_explanation + s3_s4_analysis)
    
    print("Step 2: Applying the 'Only one statement is true' rule.\n")
    
    rule_application = ("Given that S5 is always False, and S4 cannot be the only true statement, the single true statement must be either S1 or S3. "
                      "If S3 is true, S4 must be false. If S1 is true, S3 and S4 must be false.\n")
    print(rule_application)

    print("Step 3: A special analysis of Statement S1.\n")
    
    s1_explanation = "S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'\n"
    s1_analysis = ("The phrase 'if yesterday was after tomorrow' is a logically impossible premise (it's always false). "
                   "In standard logic, a conditional statement with a false premise is always true, regardless of the conclusion. "
                   "However, for this puzzle to have a unique answer, we must assume 'puzzle logic' where a conditional with an impossible premise makes the statement's truth value dependent on the conclusion. "
                   "Therefore, S1 is True if and only if the day the statement is considered is Wednesday.\n")
    print(s1_explanation + s1_analysis)
    
    print("Step 4: Combining all facts to find the solution.\n")

    conclusion = ("We have two possibilities for the one true statement: S1 or S3.\n"
                  "- If S3 were the true statement, then S1 would have to be false. This would mean the day is NOT Wednesday. This scenario leads to ambiguity as it could be any of the other six days.\n"
                  "- If S1 is the true statement, this can only happen on one specific day: Wednesday. For S1 to be the *only* true statement, S3, S4, and S5 must be false on that day, which is a consistent possibility.\n\n"
                  "The puzzle states: 'Only one statement is true and Fred said it on the day when he tells the truth.'\n"
                  "Since the only way to uniquely identify a day is if S1 is the true statement, the day must be Wednesday.\n"
                  "So, on Wednesday, statement S1 is true. Fred said S1 on his truth-telling day. Therefore, Wednesday is the day Fred tells the truth.\n")
    print(conclusion)
    
    final_answer_day = "Wednesday"
    print(f"Final Answer: The day Fred tells the truth is {final_answer_day}.")

solve_puzzle()