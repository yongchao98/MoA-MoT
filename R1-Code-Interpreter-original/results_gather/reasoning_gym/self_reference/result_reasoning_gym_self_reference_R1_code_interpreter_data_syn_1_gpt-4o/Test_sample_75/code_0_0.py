def check_statements():
    # Possible number of true statements
    for true_count in range(8):
        false_count = 7 - true_count
        
        # Check each statement
        s1 = (true_count >= 7)
        s2 = (false_count <= 3)
        s3 = (true_count == 4)
        s4 = (false_count == 7)
        s5 = (s3 != s4)
        s6 = (true_count in [2, 3, 5, 7])  # Prime numbers
        s7 = (false_count in [4, 6])  # Composite numbers
        
        # Count the number of true statements
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
        
        # Check if the number of true statements matches the true_count
        if true_statements == true_count:
            return 1  # Only one solution is possible
    
    return 0  # No solution found

print(check_statements())