def solve_freds_lie_puzzle():
    """
    Solves the logic puzzle to determine which day Fred tells the truth.
    This function analyzes the statements and simulates the conditions for each day of the week.
    """
    days_of_week = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
    
    # Let's assign numbers to days for easier calculation (Monday=0, Sunday=6)
    day_map = {day: i for i, day in enumerate(days_of_week)}
    wednesday_index = day_map["Wednesday"]

    def is_s1_true(today_index):
        """
        S1: "if yesterday was after tomorrow, it would be Wednesday."
        This is a riddle. "after tomorrow" means Today + 2. "yesterday" means Today - 1.
        The shift is (Today + 2) - (Today - 1) = +3 days.
        The statement is true if Today + 3 days = Wednesday.
        """
        # The equation is (today_index + 3) % 7 == wednesday_index
        # So, today_index == (wednesday_index - 3 + 7) % 7
        target_day_index = (wednesday_index - 3 + 7) % 7
        return today_index == target_day_index

    # We need to find the day where exactly one statement is true.
    # We determined S4 and S5 are always false.
    # The puzzle only works if S3 is also false, leaving S1 as the only possibility.
    s3_is_true = False 
    s4_is_true = False
    s5_is_true = False

    print("Analyzing the statements step-by-step:")
    print("S5: 'My best friend is older than me. He is 18 and I am 19.' -> This is a contradiction, so S5 is always False.")
    print("S4: 'The number of the males is equal to number of the females.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("Logic: If S4 is True, the total number of friends is even, so S3 must be True. This makes two true statements.")
    print("Rule: 'Only one statement is true.' -> Therefore, S4 must be False.")
    
    print("\nAnalyzing S1: 'if yesterday was after tomorrow, it would be Wednesday.'")
    shift = 3
    print(f"This is a riddle implying a time shift. The shift from 'yesterday' to 'after tomorrow' is {shift} days.")
    target_day_index = (day_map["Wednesday"] - shift + 7) % 7
    s1_truth_day = days_of_week[target_day_index]
    print(f"The equation is: Today + {shift} days = Wednesday.")
    print(f"Solving for Today: Wednesday - {shift} days = {s1_truth_day}.")
    print(f"Therefore, S1 is only true on {s1_truth_day}.\n")
    
    print("To ensure a unique solution, we deduce that S3 must be False.")
    print("If S3 were True, any day S1 is False (Mon-Sat) would be a valid day, which isn't a unique answer.")
    print("So, the only true statement can possibly be S1.\n")

    solution_day = None
    
    print("Finding the single day where the condition 'only one statement is true' holds:")
    for i, day in enumerate(days_of_week):
        s1_is_true_today = is_s1_true(i)
        
        true_statements = []
        if s1_is_true_today:
            true_statements.append("S1")
        if s3_is_true:
            true_statements.append("S3")
        if s4_is_true:
            true_statements.append("S4")
        if s5_is_true:
            true_statements.append("S5")

        # Check if the puzzle's main condition is met on this day
        if len(true_statements) == 1:
            print(f"- On {day}: S1 is {s1_is_true_today}. The set of true statements is {true_statements}. This fits the condition.")
            solution_day = day
            # Since we reasoned for a unique solution, we can stop here.
            break
        else:
            print(f"- On {day}: S1 is {s1_is_true_today}. The set of true statements is {true_statements}. This does not fit.")
            
    if solution_day:
        print(f"\nConclusion: The unique day that satisfies all puzzle conditions is {solution_day}.")
        print(f"Fred tells the truth on {solution_day}.")

solve_freds_lie_puzzle()
<<<G>>>