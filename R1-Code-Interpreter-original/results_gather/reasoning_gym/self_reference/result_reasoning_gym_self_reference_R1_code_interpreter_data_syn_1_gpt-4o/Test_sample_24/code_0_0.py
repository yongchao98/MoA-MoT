from sympy import isprime

def check_statements():
    solutions = 0
    for true_count in range(8):  # 0 to 7 true statements
        false_count = 7 - true_count
        
        # Evaluate each statement
        s1 = true_count >= 5
        s2 = false_count <= 5
        s3 = true_count == 5
        s4 = false_count == 2
        s5 = (s3 or s4) and not (s3 and s4)
        s6 = isprime(true_count)
        s7 = false_count in [4, 6]  # Composite numbers less than 7
        
        # Check if the number of true statements matches the count
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
        
        if true_statements == true_count:
            solutions += 1
    
    print(solutions)

check_statements()