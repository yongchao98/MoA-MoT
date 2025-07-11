def solve_freds_riddle():
    """
    This function solves the logic puzzle about Fred's lying day.
    It analyzes each statement to determine its truth value and then
    deduces the day based on the puzzle's rules.
    """

    # --- Step 1: Analyze the statements to find the single true one ---

    print("Step 1: Analyzing the statements to find the one true statement.")

    # Statement 5 Analysis
    s5_is_true = False
    my_age = 19
    friend_age = 18
    # The statement claims the friend is older, but 18 is not > 19.
    # It's a self-contradiction.
    print(f"S5: 'My best friend is older than me. He is {friend_age} and I am {my_age}.'")
    print(f"Analysis: This is a contradiction because {friend_age} is not older than {my_age}. So, S5 is FALSE.")
    print("-" * 20)

    # Statement 4 and 3 Analysis
    s4_is_true = False
    # If S4 ('males == females') were true, the total number of friends would be M+M = 2M, which is always even.
    # This would make S3 ('total is even') also true.
    # The puzzle states only one statement is true, so S4 cannot be true.
    print("S4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("Analysis: If S4 were true, S3 would also be true. Since only one statement can be true, S4 must be FALSE.")
    print("-" * 20)

    # Statement 1 Analysis
    s1_is_true = True
    # The statement is "if P, then Q", where P is "yesterday was after tomorrow".
    # The premise P is a logical impossibility, so it is always false.
    # In logic, an "if-then" statement with a false premise is always true.
    print("S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("Analysis: The premise 'yesterday was after tomorrow' is impossible (always false). A conditional statement with a false premise is always true. So, S1 is TRUE.")
    print("-" * 20)

    # --- Step 2: Synthesize and find the answer ---

    print("Step 2: Synthesizing the results.")
    print("We have determined that S1 is the only true statement.")
    print("The puzzle states Fred said the true statement (S1) on a day he tells the truth.")
    print("\nThe statement S1 is: '...if yesterday was after tomorrow, it would be Wednesday.'")
    print("This statement is always true, no matter what day it is.")
    print("The puzzle's trick is to connect the content of the only true statement to the answer.")
    print("The statement mentions only one day of the week: Wednesday.")
    print("\nTherefore, the day in question is Wednesday.")

    # The final equation is a representation of the logic.
    # Let's represent the truth values as 1 for True and 0 for False.
    # The puzzle is solved by finding the one true statement.
    s1_val = 1
    s3_val = 0
    s4_val = 0
    s5_val = 0
    print("\nFinal Equation of Truth Values:")
    print(f"S1({s1_val}) + S3({s3_val}) + S4({s4_val}) + S5({s5_val}) = {s1_val + s3_val + s4_val + s5_val} (The one true statement is S1)")


solve_freds_riddle()
<<<C>>>