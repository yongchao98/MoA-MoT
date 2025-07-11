def solve_puzzle():
    """
    This function explains the step-by-step solution to the logic puzzle.
    """
    print("Let's solve the puzzle step by step:")
    
    print("\nStep 1: Analyze Statement S5.")
    print("S5: 'My best friend is older than me. He is 18 and I am 19.'")
    print("This is a contradiction. If Fred is 19 and his friend is 18, his friend is younger.")
    print("Result: S5 is FALSE.")

    print("\nStep 2: Analyze the relationship between S3 and S4.")
    print("S4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("If S4 is true (males = females), then the total number of friends is males + females = 2 * males, which is always an even number.")
    print("This means that if S4 is true, S3 must also be true.")
    
    print("\nStep 3: Apply the 'Only one statement is true' rule.")
    print("The puzzle states only one statement is true. Since S4 being true would force S3 to also be true (making two true statements), S4 cannot be the true statement.")
    print("Result: S4 is FALSE.")

    print("\nStep 4: Analyze Statement S1.")
    print("S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("The condition 'if yesterday was after tomorrow' is a logical impossibility, so it is always false.")
    print("In logic, a conditional statement with a false premise is always considered TRUE.")
    print("Assuming the speaker's name is indeed Fred, the entire statement is true.")
    print("Result: S1 is TRUE.")

    print("\nStep 5: Identify the single true statement.")
    print("We have found one true statement (S1) and two false statements (S4, S5).")
    print("Since only one statement can be true, S1 must be the one. This also means S3 must be false.")
    
    print("\nStep 6: Determine the day of the week.")
    print("The puzzle states that Fred said the true statement (S1) on a day he tells the truth.")
    print("The true statement is: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("While this statement is logically true on any day, the only clue within the statement that points to a specific day is the word 'Wednesday'.")
    print("In this type of logic puzzle, such a clue is typically the key to the solution.")
    
    print("\nFinal Conclusion: The day Fred tells the truth is Wednesday.")

solve_puzzle()