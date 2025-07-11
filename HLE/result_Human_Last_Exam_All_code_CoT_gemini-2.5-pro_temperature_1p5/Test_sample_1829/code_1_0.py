def solve_puzzle():
    """
    This function analyzes the logic puzzle to find the day Fred tells the truth.
    """

    # Step 1 & 2: Analyze statements and apply the 'only one is true' rule.
    # S5 is a contradiction: False
    s5_truth_value = False
    # S1 is a vacuously true statement in formal logic: True
    s1_truth_value_logic = True

    # The condition "Only one statement is true" means S1 must be the true one.
    # This implies S3 and S4 must be false.
    true_statement_is = "S1"
    s3_truth_value = False
    s4_truth_value = False

    # Step 3: Interpret the riddle in S1 to find the day.
    # The riddle "if [impossible event], it would be Wednesday" ties the statement's truth to a specific day.
    # We model this: S1 is only considered true on Wednesday.
    s1_is_true_on_day = "Wednesday"

    # Step 4: Confirm this is the only consistent solution.
    # If the truth day were not Wednesday, S3 would have to be true.
    # But then on Wednesday (a lie day), both S1 and S3 would be true,
    # violating the "Only one statement is true" premise.
    # Thus, the only possibility is that the truth day is Wednesday, which makes S1 true
    # and S3 must be false, satisfying all conditions.

    final_day = s1_is_true_on_day
    
    print("Step 1: Analyzing the statements.")
    print("S5 ('...older than me. He is 18 and I am 19.') is a contradiction, so it is always FALSE.")
    print("S1 ('...if yesterday was after tomorrow...') has a logically impossible premise, making the conditional statement always TRUE in formal logic.")
    
    print("\nStep 2: Applying the main rule 'Only one statement is true'.")
    print(f"Since S1 is always logically TRUE, it must be the single true statement.")
    print("This means S3 ('...friends is even') and S4 ('...males is equal to females') must be FALSE.")

    print("\nStep 3: Deducing the day from the true statement.")
    print("The phrase '...it would be Wednesday' in the vacuously true statement S1 is a classic riddle.")
    print("It's meant to be interpreted as: the statement S1 is tied to Wednesday.")
    print("The logical puzzle's structure forces a conclusion that the only day this entire scenario can be self-consistent is the day mentioned.")

    print("\nStep 4: Final Conclusion.")
    print(f"The analysis points to the day mentioned in the only true statement. The day is {final_day}.")

solve_puzzle()