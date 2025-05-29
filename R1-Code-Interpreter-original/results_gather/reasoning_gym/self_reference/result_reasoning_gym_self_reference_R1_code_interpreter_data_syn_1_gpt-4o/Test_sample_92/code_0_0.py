from sympy import isprime

def check_statements():
    solutions = 0
    for T in range(8):
        F = 7 - T
        # Check each statement
        s1 = T >= 4
        s2 = F <= 1
        s3 = T == 0
        s4 = F == 3
        s5 = (s3 and not s4) or (s4 and not s3)
        s6 = isprime(T)
        s7 = F in [4, 6, 8]
        
        # Count the number of true statements
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
        
        # Check if the number of true statements equals T
        if true_statements == T:
            solutions += 1
    
    print(solutions)

check_statements()