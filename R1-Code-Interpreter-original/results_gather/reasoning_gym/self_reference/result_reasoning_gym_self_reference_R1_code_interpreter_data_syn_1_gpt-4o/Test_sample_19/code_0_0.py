from itertools import product

def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

# Generate all combinations of truth values for 7 statements
combinations = list(product([True, False], repeat=7))
valid_solutions = 0

for combo in combinations:
    s1, s2, s3, s4, s5, s6, s7 = combo
    
    # Check each statement's condition
    if (s1 == (sum(combo) >= 7) and
        s2 == (sum(combo) >= 6) and
        s3 == (sum(combo) == 1) and
        s4 == (sum(combo) == 0) and
        s5 == (s3 != s4) and
        s6 == is_prime(sum(combo)) and
        s7 == is_composite(7 - sum(combo))):
        valid_solutions += 1

print(valid_solutions)