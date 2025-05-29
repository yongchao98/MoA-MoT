def check_statements(truth_values):
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # Check each statement
    statements = [
        # Statement 1: At least 1 of these 7 statements are true
        true_count >= 1,
        
        # Statement 2: At most 3 of these 7 statements are false
        false_count <= 3,
        
        # Statement 3: Exactly 4 of these 7 statements are true
        true_count == 4,
        
        # Statement 4: Exactly 3 of these 7 statements are false
        false_count == 3,
        
        # Statement 5: Either Statement 3 or Statement 4 is true, but not both
        truth_values[2] != truth_values[3] and (truth_values[2] or truth_values[3]),
        
        # Statement 6: The number of true statements is a prime number
        true_count in [2, 3, 5, 7],
        
        # Statement 7: The number of false statements is a composite number
        false_count in [4, 6]  # composite numbers <= 7
    ]
    
    # Check if the truth values match the actual evaluation of statements
    return all(truth_values[i] == statements[i] for i in range(7))

# Generate all possible combinations and check each
valid_solutions = []
for i in range(2**7):
    # Convert number to binary representation of 7 bits
    truth_values = [bool((i >> j) & 1) for j in range(7)]
    if check_statements(truth_values):
        valid_solutions.append(truth_values)

print(f"Number of valid solutions: {len(valid_solutions)}")
if valid_solutions:
    print("\nValid solutions:")
    for solution in valid_solutions:
        print([int(x) for x in solution])