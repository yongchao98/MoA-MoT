def check_solution(true_count):
    # Given true_count statements are true, check if it's consistent
    false_count = 7 - true_count
    
    # Initialize all statements as false
    statements = [False] * 7
    
    # Statement 1: At least 7 are true
    statements[0] = (true_count >= 7)
    
    # Statement 2: At most 7 are false
    statements[1] = (false_count <= 7)
    
    # Statement 3: Exactly 5 are true
    statements[2] = (true_count == 5)
    
    # Statement 4: Exactly 2 are false
    statements[3] = (false_count == 2)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    statements[4] = (statements[2] != statements[3] and (statements[2] or statements[3]))
    
    # Statement 6: Number of true statements is prime
    statements[5] = true_count in {2, 3, 5, 7}
    
    # Statement 7: Number of false statements is composite
    statements[6] = false_count in {4, 6}
    
    # Count how many statements are actually true
    actual_true = sum(statements)
    
    # Solution is valid if actual_true equals assumed true_count
    return actual_true == true_count

# Check all possible numbers of true statements (0 to 7)
solutions = []
for true_count in range(8):
    if check_solution(true_count):
        solutions.append(true_count)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for true_count in solutions:
        print(f"Number of true statements: {true_count}")
        print(f"Number of false statements: {7-true_count}")