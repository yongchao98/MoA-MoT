def solve_puzzle():
    """
    This function solves the logic puzzle by analyzing the statements step by step.
    """
    s1_text = "S1) My name is Fred; and if yesterday was after tomorrow, it would be Wednesday."
    s3_text = "S3) The total number of my friend in my class is even."
    s4_text = "S4) The number of the males is equal to number of the females in my class."
    s5_text = "S5) My best friend is older than me. He is 18 and I am 19."

    print("Step 1: Analyze the inherent truth of each statement.")
    
    # Analyzing S5
    friend_age = 18
    my_age = 19
    is_s5_true = friend_age > my_age
    print(f"Statement S5 claims that the friend ({friend_age}) is older than me ({my_age}).")
    print(f"The check is: {friend_age} > {my_age}, which is {is_s5_true}.")
    print("Therefore, S5 is definitively False.\n")

    # Analyzing S4 and S3
    print("Statement S4 says the number of males (M) equals the number of females (F).")
    print("Statement S3 says the total number of friends (M + F) is even.")
    print("If S4 is true, then M = F. The total is M + F = 2 * M, which is always even.")
    print("This means if S4 is true, S3 must also be true.")
    print("However, the puzzle states *only one* statement can be true.")
    print("Therefore, S4 cannot be the single true statement, so S4 must be False.\n")

    # Analyzing S1
    print("Step 2: Identify the puzzle's trick in statement S1.")
    print(f"Statement S1 is: '{s1_text}'")
    print("Logically, 'if yesterday was after tomorrow' is an impossible premise, making the whole conditional 'vacuously true'.")
    print("If S1 is always true, it must be the single true statement. This doesn't help us find the day.")
    print("The puzzle's trick is likely that a statement 'if [absurdity], then Q' is a convoluted way of asserting 'NOT Q'.")
    print("So, we interpret S1 as meaning: 'My name is Fred' AND 'Today is NOT Wednesday'.\n")

    # Finding the day
    print("Step 3: Find the only day that produces a valid scenario.")
    print("We need to find a day where only one statement is true.")
    print("The candidates for the single true statement are S1 or S3.")
    print("- If S1 ('Today is NOT Wednesday') were the true statement, the day could be any of 6 days. This is not a unique solution.")
    print("- Let's test if S3 can be the single true statement.")
    
    # This means S1 must be false for S3 to be the *only* true statement.
    print("  For S3 to be the ONE true statement, S1 MUST be false.")
    s1_is_false_condition = "the statement 'Today is NOT Wednesday' must be false"
    print(f"  For S1 to be false, {s1_is_false_condition}.")
    final_conclusion_day = "Wednesday"
    print(f"  This is only possible if today is {final_conclusion_day}.")
    
    print("\nStep 4: Final Conclusion")
    print("Let's assume today is Wednesday:")
    print(f"- S1 ('Today is not Wednesday') becomes False.")
    print(f"- S3 can now be the single True statement.")
    print(f"- S4 is False.")
    print(f"- S5 is False.")
    print("This scenario (S3 is True, others are False) works only on one day.")
    print(f"The day Fred tells the truth is {final_conclusion_day}.")

solve_puzzle()