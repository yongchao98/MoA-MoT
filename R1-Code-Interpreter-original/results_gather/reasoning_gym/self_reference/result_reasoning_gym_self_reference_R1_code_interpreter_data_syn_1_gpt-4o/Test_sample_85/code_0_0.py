from itertools import combinations

def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(true_statements):
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    # Check each statement
    s1 = num_true >= 4
    s2 = num_false <= 5
    s3 = num_true == 4
    s4 = num_false == 7
    s5 = (s3 or s4) and not (s3 and s4)
    s6 = is_prime(num_true)
    s7 = is_composite(num_false)
    
    # Determine which statements are true
    true_set = {1, 2, 3, 5}
    false_set = {4, 6, 7}
    
    return (s1, s2, s3, s4, s5, s6, s7) == (s1 in true_set, s2 in true_set, s3 in true_set, s4 in true_set, s5 in true_set, s6 in true_set, s7 in true_set)

# Check all combinations of true statements
possible_solutions = 0
for true_statements in combinations(range(1, 8), 4):
    if check_statements(true_statements):
        possible_solutions += 1

print(possible_solutions)