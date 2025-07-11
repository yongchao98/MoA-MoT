def solve_fred_logic_puzzle():
    """
    This function solves the logic puzzle by analyzing the truth value of each statement
    based on the given rules.
    """

    print("Step 1: Analyzing the inherent truth value of each statement.")
    
    # Statement S1: "My name is Fred; and if yesterday was after tomorrow, it would be Wednesday."
    # The condition "yesterday was after tomorrow" is a logical impossibility (always False).
    # In logic, an "if P then Q" implication is considered TRUE whenever the premise P is False.
    # This is known as being vacuously true. Therefore, S1 is a logically TRUE statement.
    s1_is_true = True
    print("- Statement S1 ('...if yesterday was after tomorrow, it would be Wednesday') is a vacuously true statement. Its premise is impossible, which makes the entire logical implication TRUE.")

    # Statement S5: "My best friend is older than me. He is 18 and I am 19."
    # This statement contains an internal contradiction. If Fred is 19 and his friend is 18,
    # his friend is younger, not older. Therefore, S5 is FALSE.
    s5_is_true = False
    print("- Statement S5 ('...He is 18 and I am 19') is self-contradictory. A 19-year-old is older than an 18-year-old. The statement is FALSE.")
    
    print("\nStep 2: Analyzing the relationship between statements S3 and S4.")
    # Statement S4: "The number of the males is equal to number of the females in my class."
    # Statement S3: "The total number of my friend in my class is even."
    # If S4 is true (males = females), then the total number of friends is males + females = 2 * males.
    # Any number multiplied by 2 is even. So, if S4 is true, S3 must also be true.
    print("- If S4 ('males equals females') is true, then the total number of students must be even. Therefore, if S4 is TRUE, S3 must also be TRUE.")

    print("\nStep 3: Applying the puzzle's main rule: 'Only one statement is true'.")
    
    # We already know S1 is TRUE.
    print("- From Step 1, we know S1 is TRUE.")
    
    # If S4 were true, S3 would also be true, violating the "only one is true" rule.
    # Therefore, S4 must be false.
    s4_is_true = False
    print("- If S4 were the true statement, S3 would also be true. This would mean two statements are true, which violates the rule. Thus, S4 must be FALSE.")

    # Since S1 is the one true statement, all other statements must be false.
    s3_is_true = False
    print("- Because S1 is true and it's the ONLY true statement, S3 must also be FALSE.")
    
    print("\nSummary of Truth Values:")
    print(f"- S1 is True: {s1_is_true}")
    print(f"- S3 is False: {s3_is_true}")
    print(f"- S4 is False: {s4_is_true}")
    print(f"- S5 is False: {s5_is_true}")

    print("\nStep 4: Determining the truth-telling day.")
    # The single true statement is S1. Fred said this statement on his truth-telling day.
    # The statement itself contains the name of the day.
    final_day = "Wednesday"
    print(f"- The only true statement is S1, which Fred spoke on his truth-telling day.")
    print(f"- Statement S1 explicitly mentions the day '{final_day}'.")
    print(f"- Therefore, the day Fred tells the truth is {final_day}.")

solve_fred_logic_puzzle()