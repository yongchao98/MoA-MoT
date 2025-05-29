def check_statements_strict(truth_values):
    # Count true statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    # Check each statement
    statements = [
        true_count >= 6,                    # Statement 1
        false_count <= 5,                   # Statement 2
        true_count == 2,                    # Statement 3
        false_count == 2,                   # Statement 4
        (truth_values[2] != truth_values[3] and 
         (truth_values[2] or truth_values[3])),  # Statement 5
        is_prime(true_count),               # Statement 6
        is_composite(false_count)           # Statement 7
    ]
    
    # First check if the statements match their truth values
    if not all(truth_values[i] == statements[i] for i in range(7)):
        return False
        
    # Additional validation
    if true_count + false_count != 7:
        return False
    
    return True

valid_solutions = []
for i in range(2**7):
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements_strict(truth_values):
        valid_solutions.append(truth_values)

print(f"Number of valid solutions: {len(valid_solutions)}")
for i, solution in enumerate(valid_solutions, 1):
    print(f"\nValid Solution {i}:")
    print(f"Truth values: {solution}")
    print(f"True statements: {sum(solution)}")
    print(f"False statements: {7 - sum(solution)}")