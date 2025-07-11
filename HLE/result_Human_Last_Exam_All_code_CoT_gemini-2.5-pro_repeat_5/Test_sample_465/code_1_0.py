def solve_task():
    """
    Analyzes the statements about chromatic and orbital chromatic roots and prints the sorted string of true statements.
    """
    # Statement A: True. Real orbital chromatic roots are bounded by the greatest real chromatic root.
    # Statement B: True. Chromatic roots are not always real.
    # Statement C: True. Real chromatic roots can be negative.
    # Statement D: True. Real chromatic roots can be non-integers.
    # Statement E: False. There are no chromatic roots in the interval (0, 1).

    true_statements = ['A', 'B', 'C', 'D']
    
    # Sort the list of true statements alphabetically
    true_statements.sort()
    
    # Join the characters to form the final answer string
    answer = "".join(true_statements)
    
    print(answer)

solve_task()