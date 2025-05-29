def check_statements(true_statements):
    # Count of true and false statements
    true_count = len(true_statements)
    false_count = 7 - true_count
    
    # Helper function to check if a statement is true
    def is_true(stmt_num):
        return stmt_num in true_statements
    
    # Check each statement's validity
    valid = True
    
    # Statement 1: At least 2 are true
    if is_true(1):
        valid &= (true_count >= 2)
    else:
        valid &= (true_count < 2)
    
    # Statement 2: At most 4 are false
    if is_true(2):
        valid &= (false_count <= 4)
    else:
        valid &= (false_count > 4)
    
    # Statement 3: Exactly 4 are true
    if is_true(3):
        valid &= (true_count == 4)
    else:
        valid &= (true_count != 4)
    
    # Statement 4: Exactly 3 are false
    if is_true(4):
        valid &= (false_count == 3)
    else:
        valid &= (false_count != 3)
    
    # Statement 5: Either Statement 3 or 4 is true, but not both
    if is_true(5):
        valid &= (is_true(3) != is_true(4))
    else:
        valid &= (is_true(3) == is_true(4))
    
    # Statement 6: Number of true statements is prime
    primes = {2, 3, 5, 7}
    if is_true(6):
        valid &= (true_count in primes)
    else:
        valid &= (true_count not in primes)
    
    # Statement 7: Number of false statements is composite
    composites = {4, 6}  # only possible composites for 7 statements
    if is_true(7):
        valid &= (false_count in composites)
    else:
        valid &= (false_count not in composites)
    
    return valid

# Generate all possible combinations and show valid solutions
valid_solutions = 0
for i in range(128):  # 2^7 possibilities
    true_statements = set()
    for j in range(7):
        if (i >> j) & 1:
            true_statements.add(j + 1)
    if check_statements(true_statements):
        valid_solutions += 1
        print(f"Solution {valid_solutions}: True statements: {sorted(true_statements)}, "
              f"False statements: {sorted(set(range(1, 8)) - true_statements)}")

print(f"\nTotal number of valid solutions: {valid_solutions}")