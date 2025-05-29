def check_statements():
    possible_solutions = 0
    for true_statements in range(8):
        false_statements = 7 - true_statements
        # Check each statement
        s1 = true_statements >= 6
        s2 = false_statements <= 5
        s3 = true_statements == 2
        s4 = false_statements == 2
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = true_statements in [2, 3, 5, 7]
        s7 = false_statements in [4, 6]
        
        # Check if all conditions are satisfied
        if s3 and s5 and s6 and s7:
            possible_solutions += 1
    return possible_solutions

print(check_statements())