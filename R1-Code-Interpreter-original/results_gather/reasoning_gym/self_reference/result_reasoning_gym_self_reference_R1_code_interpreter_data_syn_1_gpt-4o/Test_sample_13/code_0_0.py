def check_statements():
    # Possible number of true statements
    for true_count in range(8):
        false_count = 7 - true_count
        
        # Check each statement
        s1 = (true_count >= 7)
        s2 = (false_count <= 2)
        s3 = (true_count == 0)
        s4 = (false_count == 4)
        s5 = (s3 != s4)
        s6 = (true_count in [2, 3, 5, 7])
        s7 = (false_count in [4, 6])
        
        # Count how many statements are true
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
        
        # Check if the number of true statements matches the true_count
        if true_statements == true_count:
            print(f"True statements: {true_count}, False statements: {false_count}")

check_statements()