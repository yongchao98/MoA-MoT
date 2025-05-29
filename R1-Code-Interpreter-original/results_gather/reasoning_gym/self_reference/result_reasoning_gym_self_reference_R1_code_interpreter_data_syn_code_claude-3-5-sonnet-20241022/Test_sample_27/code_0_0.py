def check_statements(combination):
    # Count true and false statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check each statement's consistency
    consistent = [False] * 7
    
    # Statement 1: "At least 2 of these 7 statements are true"
    consistent[0] = (true_count >= 2) == combination[0]
    
    # Statement 2: "At most 0 of these 7 statements are false"
    consistent[1] = (false_count <= 0) == combination[1]
    
    # Statement 3: "Exactly 1 of these 7 statements are true"
    consistent[2] = (true_count == 1) == combination[2]
    
    # Statement 4: "Exactly 7 of these 7 statements are false"
    consistent[3] = (false_count == 7) == combination[3]
    
    # Statement 5: "Either Statement 3 or Statement 4 is true, but not both"
    consistent[4] = ((combination[2] != combination[3]) and (combination[2] or combination[3])) == combination[4]
    
    # Statement 6: "The number of true statements is a prime number"
    primes = {2, 3, 5, 7}
    consistent[5] = (true_count in primes) == combination[5]
    
    # Statement 7: "The number of false statements is a composite number"
    composites = {4, 6}  # only possible composite numbers for 7 statements
    consistent[6] = (false_count in composites) == combination[6]
    
    return all(consistent)

# Generate all possible combinations and check them
solutions = []
for i in range(2**7):
    combination = [(i >> j) & 1 for j in range(7)]
    if check_statements(combination):
        solutions.append(combination)

print(f"Number of solutions: {len(solutions)}")
if solutions:
    print("\nSolutions found:")
    for sol in solutions:
        print(f"Truth values: {sol}")
        print(f"Number of true statements: {sum(sol)}")
        print(f"Number of false statements: {7-sum(sol)}")
        print()