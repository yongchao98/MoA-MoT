def solve_puzzle():
    """
    This function solves the logic puzzle by analyzing each statement
    and applying the given constraints.
    """
    
    # Introduction to the analysis
    print("Analyzing the puzzle step by step:")
    print("The core rule is: Exactly one of the four statements is true.")
    print("-" * 30)

    # Step 1: Analyze Statement S5
    # S5: "My best friend is older than me. He is 18 and I am 19."
    fred_age = 19
    friend_age = 18
    is_friend_older = friend_age > fred_age
    print("Step 1: Analyze Statement S5.")
    print(f"S5 claims a friend of age {friend_age} is older than Fred, age {fred_age}.")
    print(f"The comparison is: {friend_age} > {fred_age}, which is {is_friend_older}.")
    print("Conclusion: Statement S5 is definitively False.\n")

    # Step 2: Analyze the relationship between S4 and S3
    # S4: "The number of the males is equal to number of the females in my class."
    # S3: "The total number of my friend in my class is even."
    print("Step 2: Analyze the relationship between S3 and S4.")
    print("S4 says: males = females.")
    print("S3 says: total_friends is even.")
    print("If S4 is true, then total_friends = males + females = 2 * males.")
    print("A number multiplied by 2 is always even. Therefore, if S4 is true, S3 must also be true.")
    print("-" * 30)

    # Step 3: Apply the "only one is true" rule to S3 and S4
    print("Step 3: Apply the 'only one true' rule.")
    print("If S4 were true, S3 would also be true. This would mean two statements are true, which is not allowed.")
    print("Conclusion: Statement S4 must be False.\n")
    
    # Step 4: Analyze Statement S1
    # S1: "My name is Fred; and if yesterday was after tomorrow, it would be Wednesday."
    print("Step 4: Analyze Statement S1.")
    print("The statement has two parts: 'My name is Fred' (A) and 'if yesterday was after tomorrow, it would be Wednesday' (B).")
    print("Part A is a given premise, so it's True.")
    print("Part B is a conditional statement. The premise 'yesterday was after tomorrow' is a logical impossibility, making it False.")
    print("In formal logic, any conditional statement with a False premise ('if False, then ...') is always considered True.")
    print("Therefore, Part B is True.")
    print("The full statement S1 is 'True and True', which makes the entire statement True.")
    print("Conclusion: Statement S1 is logically True.\n")

    # Step 5: Final summary and answer
    print("Step 5: Final determination.")
    print("We have concluded:")
    print("- S1 is True.")
    print("- S5 is False.")
    print("- S4 is False.")
    print("Since only one statement can be true, S3 must also be False to satisfy the rule.")
    print("The single true statement is S1.\n")
    
    print("Step 6: The Answer.")
    print("The only true statement is S1, which states: '...if yesterday was after tomorrow, it would be Wednesday.'")
    print("The question asks which day Fred tells the truth. The single true statement explicitly mentions 'Wednesday'.")
    print("This is the answer pointed to by the logic of the puzzle.")

solve_puzzle()
<<<C>>>