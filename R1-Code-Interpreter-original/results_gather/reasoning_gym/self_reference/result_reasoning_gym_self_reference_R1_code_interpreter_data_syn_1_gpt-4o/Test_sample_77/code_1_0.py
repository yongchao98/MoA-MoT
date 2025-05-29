from sympy import isprime

def check_statements():
    possible_solutions = 0
    for true_count in range(8):
        false_count = 7 - true_count
        # Check each statement
        s1 = (true_count >= 7)
        s2 = (false_count <= 7)
        s3 = (true_count == 5)
        s4 = (false_count == 0)
        s5 = (s3 != s4)
        s6 = isprime(true_count)
        s7 = (false_count in [4, 6])
        
        # Count valid solutions
        if s5 and s6 and s7:
            possible_solutions += 1
    return possible_solutions

print(check_statements())