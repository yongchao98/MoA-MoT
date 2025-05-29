def check_solution(pattern):
    # Count true and false statements based on the current pattern
    true_count = sum(1 for x in pattern)
    false_count = 7 - true_count
    
    # Evaluate what each statement claims
    statement_claims = [False] * 7
    
    # Statement 1: At least 4 are true
    statement_claims[0] = (true_count >= 4)
    
    # Statement 2: At most 2 are false
    statement_claims[1] = (false_count <= 2)
    
    # Statement 3: Exactly 1 is true
    statement_claims[2] = (true_count == 1)
    
    # Statement 4: All 7 are false
    statement_claims[3] = (true_count == 0)
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    # This means exactly one of statements 3 and 4 must be true in the pattern
    statement_claims[4] = (pattern[2] ^ pattern[3])  # XOR operation
    
    # Statement 6: Number of true statements is prime
    statement_claims[5] = true_count in {2, 3, 5, 7}
    
    # Statement 7: Number of false statements is composite
    statement_claims[6] = false_count in {4, 6}
    
    # Check if the actual truth values match what each statement claims
    matches = [pattern[i] == statement_claims[i] for i in range(7)]
    
    # For debugging, print detailed information
    if sum(matches) >= 5:  # Only print near-matches for debugging
        print(f"\nPattern: {[int(x) for x in pattern]}")
        print(f"Claims:  {[int(x) for x in statement_claims]}")
        print(f"Matches: {[int(x) for x in matches]}")
        print(f"True count: {true_count}, False count: {false_count}")
    
    return all(matches)

# Try all possible combinations
solutions = []
for i in range(2**7):
    # Convert number to binary pattern
    pattern = [bool((i >> j) & 1) for j in range(7)]
    if check_solution(pattern):
        solutions.append(pattern)

print(f"\nTotal number of solutions found: {len(solutions)}")
if solutions:
    print("\nValid solutions:")
    for sol in solutions:
        print([int(x) for x in sol])