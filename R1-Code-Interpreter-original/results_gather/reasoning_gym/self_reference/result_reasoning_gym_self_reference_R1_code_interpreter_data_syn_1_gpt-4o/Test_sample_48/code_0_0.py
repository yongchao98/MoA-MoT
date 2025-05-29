from itertools import product

def is_prime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

# Generate all combinations of truth values for 7 statements
combinations = product([True, False], repeat=7)
valid_solutions = 0

for combo in combinations:
    s1, s2, s3, s4, s5, s6, s7 = combo
    true_count = sum(combo)
    false_count = 7 - true_count
    
    # Check each statement's condition
    if (s1 == (true_count >= 7) and
        s2 == (false_count <= 2) and
        s3 == (true_count == 2) and
        s4 == (false_count == 4) and
        s5 == (s3 != s4) and
        s6 == is_prime(true_count) and
        s7 == is_composite(false_count)):
        valid_solutions += 1

print(valid_solutions)