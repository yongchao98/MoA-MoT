def check_statements(truth_values):
    # Count true and false statements
    true_count = sum(truth_values)
    false_count = 7 - true_count
    
    # First, calculate what each statement claims
    claim1 = true_count >= 4
    claim2 = false_count <= 5
    claim3 = true_count == 0
    claim4 = false_count == 3
    claim5 = bool(truth_values[2]) != bool(truth_values[3])  # XOR of statements 3 and 4
    claim6 = true_count in {2, 3, 5, 7}  # prime numbers
    claim7 = false_count in {4, 6}  # composite numbers possible with 7 statements
    
    # Then check if each truth value matches its claim
    consistent = (
        truth_values[0] == claim1 and
        truth_values[1] == claim2 and
        truth_values[2] == claim3 and
        truth_values[3] == claim4 and
        truth_values[4] == claim5 and
        truth_values[5] == claim6 and
        truth_values[6] == claim7
    )
    
    if consistent:
        print(f"Found solution:")
        print(f"Truth values: {truth_values}")
        print(f"True count: {true_count}")
        print(f"Claims: {[claim1, claim2, claim3, claim4, claim5, claim6, claim7]}")
        print("---")
    
    return consistent

# Test all possible combinations
solutions = []
for i in range(2**7):
    truth_values = [(i >> j) & 1 == 1 for j in range(7)]
    if check_statements(truth_values):
        solutions.append(truth_values)

print(f"Total number of solutions: {len(solutions)}")