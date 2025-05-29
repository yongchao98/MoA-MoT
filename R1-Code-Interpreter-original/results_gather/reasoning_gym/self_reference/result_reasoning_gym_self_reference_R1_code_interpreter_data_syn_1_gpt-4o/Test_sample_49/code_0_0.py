from sympy import isprime

def check_statements():
    solutions = 0
    for true_count in range(8):
        false_count = 7 - true_count
        # Check each statement
        s1 = (true_count >= 7)
        s2 = (false_count <= 2)
        s3 = (true_count == 5)
        s4 = (false_count == 2)
        s5 = (s3 != s4)
        s6 = isprime(true_count)
        s7 = (false_count in [4, 6])
        
        # Count the number of true statements
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
        
        # Check if the number of true statements matches the true_count
        if true_statements == true_count:
            solutions += 1
    
    return solutions

print(check_statements())