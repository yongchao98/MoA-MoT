from itertools import product

def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_statements(truth_values):
    # Unpack truth values
    S1, S2, S3, S4, S5, S6, S7 = truth_values
    
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # Check each statement
    if S1 and true_count < 7:
        return False
    if S2 and false_count > 5:
        return False
    if S3 and true_count != 4:
        return False
    if S4 and false_count != 3:
        return False
    if S5 and not ((S3 and not S4) or (S4 and not S3)):
        return False
    if S6 and not is_prime(true_count):
        return False
    if S7 and not is_composite(false_count):
        return False
    
    return True

# Generate all possible truth value combinations for 7 statements
possible_solutions = [truth_values for truth_values in product([True, False], repeat=7) if check_statements(truth_values)]

# Output the number of possible solutions
print(len(possible_solutions))