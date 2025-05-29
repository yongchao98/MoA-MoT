def check_statements():
    # Number of statements
    total_statements = 7
    
    # Possible number of true statements
    possible_true_counts = [5, 6, 7]
    
    # Check each possible number of true statements
    for T in possible_true_counts:
        F = total_statements - T
        
        # Check each statement
        statement_1 = (T >= 5)
        statement_2 = (F <= 3)
        statement_3 = (T == 5)
        statement_4 = (F == 5)
        statement_5 = (statement_3 != statement_4)
        statement_6 = (T in [2, 3, 5, 7])  # Prime numbers
        statement_7 = (F in [4, 6])  # Composite numbers
        
        # Check if all conditions are satisfied
        if (statement_1 and statement_2 and statement_3 and not statement_4 and
            statement_5 and statement_6 and not statement_7):
            return 1  # Only one solution is possible
    
    return 0  # No solution found

# Output the number of possible solutions
print(check_statements())