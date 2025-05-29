def check_statements():
    # Possible number of true statements
    for true_count in range(8):
        false_count = 7 - true_count
        
        # Check each statement
        s1 = true_count >= 4
        s2 = false_count <= 2
        s3 = true_count == 5
        s4 = false_count == 6
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = true_count in [2, 3, 5, 7]  # Prime numbers
        s7 = false_count in [4, 6]  # Composite numbers
        
        # Count the number of true statements
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
        
        # Check if the configuration is consistent
        if true_statements == true_count:
            print(f"Consistent solution found with {true_count} true statements and {false_count} false statements.")

check_statements()