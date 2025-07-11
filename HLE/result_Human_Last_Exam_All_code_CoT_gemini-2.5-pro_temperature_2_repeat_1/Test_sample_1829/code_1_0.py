def solve_puzzle():
    """
    This function solves the logic puzzle by evaluating each statement
    to find the single true one and then determines the day of the week.
    """
    
    # Step 1: Analyze the statements to find the single true statement.
    print("Step 1: Analyzing the statements to find the single true statement.")

    # Step 2: Evaluate S5
    s5_is_older = 18 > 19
    print("\nStep 2: Evaluating S5: 'My best friend is older than me. He is 18 and I am 19.'")
    print(f"This implies 18 > 19, which is {s5_is_older}.")
    print("Statement S5 is a self-contradiction, so it must be False.")

    # Step 3: Evaluate S4 and S3
    print("\nStep 3: Evaluating S4 and S3.")
    print("S4 states the number of males equals the number of females.")
    print("If S4 were true, the total number of friends (males + females) would be even.")
    print("This would make S3 ('The total number of my friend in my class is even') also true.")
    print("This would violate the rule that 'Only one statement is true'.")
    print("Therefore, statement S4 must be False.")

    # Step 4: Evaluate S1
    print("\nStep 4: Evaluating S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    premise_is_false = True # "yesterday was after tomorrow" is an impossible condition.
    print("The condition 'if yesterday was after tomorrow...' is impossible, so its premise is False.")
    # In logic, an implication (if P then Q) with a false premise (P) is always true.
    implication_is_true = True
    print(f"A logical implication with a False premise is always True.")
    my_name_is_fred = True
    s1_is_true = my_name_is_fred and implication_is_true
    print(f"The statement is 'My name is Fred' (True) AND the implication (True), making S1 {s1_is_true}.")

    # Step 5: Identify the one true statement
    print("\nStep 5: Identifying the single true statement.")
    print("We have determined S1 is True, while S4 and S5 are False.")
    print("Because only one statement can be true, S1 must be the one. This forces S3 to be False.")

    # Step 6: Determine the day
    print("\nStep 6: Determining the day of the week.")
    true_statement = "S1"
    day_mentioned_in_s1 = "Wednesday"
    print(f"The one true statement is {true_statement}, which mentions the day '{day_mentioned_in_s1}'.")
    print("The puzzle states Fred made this statement on the day he tells the truth.")
    print("The most logical conclusion is that the day mentioned in the only true statement is the answer.")
    print(f"Therefore, Fred tells the truth on {day_mentioned_in_s1}.")

solve_puzzle()
<<<C>>>