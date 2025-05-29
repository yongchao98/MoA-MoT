def check_statements():
    # Possible number of true statements
    possible_true_counts = [2, 6, 7]
    
    # Check each possible number of true statements
    for true_count in possible_true_counts:
        false_count = 7 - true_count
        
        # Evaluate each statement based on the current true_count
        s1 = true_count >= 6
        s2 = false_count <= 1
        s3 = true_count == 6
        s4 = false_count == 5
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = true_count in [2, 3, 5, 7]
        s7 = false_count in [4, 6]
        
        # Count the number of true statements
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
        
        # Check if the number of true statements matches the true_count
        if true_statements == true_count:
            return 1  # Return 1 as a possible solution

    return 0  # Return 0 if no solution is found

# Execute the function and print the result
print(check_statements())