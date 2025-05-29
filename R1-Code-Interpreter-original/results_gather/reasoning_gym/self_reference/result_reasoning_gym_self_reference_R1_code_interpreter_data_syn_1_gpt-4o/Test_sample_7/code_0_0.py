def check_statements():
    possible_solutions = 0
    for true_count in range(8):
        false_count = 7 - true_count
        
        # Check each statement
        s1 = (true_count >= 6)
        s2 = (false_count <= 0)
        s3 = (true_count == 2)
        s4 = (false_count == 1)
        s5 = (s3 != s4)
        s6 = (true_count in [2, 3, 5, 7])
        s7 = (false_count in [4, 6])
        
        # Check if the configuration is consistent
        if (s1 + s2 + s3 + s4 + s5 + s6 + s7) == true_count:
            possible_solutions += 1
    
    return possible_solutions

print(check_statements())