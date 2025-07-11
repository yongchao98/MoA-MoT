def solve_puzzle():
    """
    This function explains the step-by-step solution to the logic puzzle.
    """
    
    print("Step 1: Analyze the statements to find the single true statement.")
    print("The puzzle states only one of S1, S3, S4, S5 is true.")
    print("-" * 20)

    print("Analyzing S5: 'My best friend is older than me. He is 18 and I am 19.'")
    friend_age = 18
    my_age = 19
    print(f"If the friend is {friend_age} and I am {my_age}, the friend is younger, not older. So, S5 is False.")
    print("-" * 20)

    print("Analyzing S4 and S3: 'The number of the males is equal to number of the females in my class' (S4) and 'The total number of my friend in my class is even' (S3).")
    print("If S4 is true, the total number of students is (males + males) = 2 * males, which is an even number.")
    print("This would mean S3 is also true.")
    print("Since only one statement can be true, S4 must be False to avoid making S3 true as well.")
    print("-" * 20)
    
    print("Analyzing S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("The premise 'yesterday was after tomorrow' is a logical impossibility, so it is always False.")
    print("In logic, an 'if P, then Q' statement is always True if the premise P is False.")
    print("Therefore, the conditional part of S1 is 'vacuously true'.")
    print("Since the first part ('My name is Fred') is also true, the entire statement S1 is True.")
    print("-" * 20)
    
    print("Conclusion of Step 1: S1 is the only true statement. S3, S4, and S5 are false.")
    print("\nStep 2: Determine the day of the week.")
    print("Fred said the one true statement (S1) on the day he tells the truth.")
    print("The statement is: 'if yesterday was after tomorrow, it would be Wednesday.'")
    print("Although this statement is logically true on any day, a truthful person would choose the conclusion 'it would be Wednesday' to be a fact.")
    print("We can therefore infer that he made the statement on the day that makes the conclusion true.")
    
    final_day = "Wednesday"
    print(f"\nFinal Answer: The conclusion of the statement is that 'it would be {final_day}'. Therefore, Fred tells the truth on {final_day}.")

solve_puzzle()
<<<C>>>