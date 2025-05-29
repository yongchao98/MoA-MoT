def check_solution(pattern):
    # Count true and false statements
    true_count = sum(1 for x in pattern if x)
    false_count = 7 - true_count
    
    # Evaluate what each statement claims to be true
    claims = [
        true_count >= 4,                    # Statement 1
        false_count <= 2,                   # Statement 2
        true_count == 1,                    # Statement 3
        false_count == 7,                   # Statement 4
        (pattern[2] ^ pattern[3]),          # Statement 5
        true_count in {2, 3, 5, 7},         # Statement 6
        false_count in {4, 6}               # Statement 7
    ]
    
    # A solution is valid if the truth value of each statement matches its claim
    matches = [claims[i] == pattern[i] for i in range(7)]
    match_count = sum(matches)
    
    if match_count == 7:  # Only print perfect matches
        print(f"\nFound valid solution:")
        print(f"Pattern: {[int(x) for x in pattern]}")
        print(f"Claims:  {[int(x) for x in claims]}")
        print(f"True count: {true_count}")
        print(f"False count: {false_count}")
        print(f"Matches: {[int(x) for x in matches]}")
        return True
    return False

# Try all possible combinations
solutions = []
for i in range(2**7):
    pattern = [bool((i >> j) & 1) for j in range(7)]
    if check_solution(pattern):
        solutions.append(pattern)

print(f"\nTotal solutions found: {len(solutions)}")