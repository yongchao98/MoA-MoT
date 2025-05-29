def check_solution(pattern):
    # First evaluate what each statement claims
    true_count = sum(1 for x in pattern if x)  # Count of True values in pattern
    false_count = 7 - true_count
    
    # What each statement claims
    claims = [
        true_count >= 4,                    # Statement 1
        false_count <= 2,                   # Statement 2
        true_count == 1,                    # Statement 3
        false_count == 7,                   # Statement 4
        (pattern[2] ^ pattern[3]),          # Statement 5 (XOR)
        true_count in {2, 3, 5, 7},         # Statement 6
        false_count in {4, 6}               # Statement 7
    ]
    
    # For a valid solution:
    # 1. Each claim must match its corresponding truth value in pattern
    # 2. The number of true statements must equal true_count
    
    # Check if claims match pattern
    valid = all(claims[i] == pattern[i] for i in range(7))
    
    # If it's close to valid, print for debugging
    if sum(claims[i] == pattern[i] for i in range(7)) >= 5:
        print(f"\nPattern: {[int(x) for x in pattern]}")
        print(f"Claims:  {[int(x) for x in claims]}")
        print(f"True count: {true_count}")
        print(f"False count: {false_count}")
        matches = [claims[i] == pattern[i] for i in range(7)]
        print(f"Matches: {[int(x) for x in matches]}")
    
    return valid

# Try all possible combinations
solutions = []
for i in range(2**7):
    pattern = [bool((i >> j) & 1) for j in range(7)]
    if check_solution(pattern):
        solutions.append(pattern)

print(f"\nTotal solutions found: {len(solutions)}")
if solutions:
    print("\nValid solutions:")
    for sol in solutions:
        print([int(x) for x in sol])
        print(f"True statements: {sum(sol)}")
        print(f"False statements: {7-sum(sol)}")