def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_composite(n):
    return n > 1 and not is_prime(n)

def check_consistency(true_count):
    # Given the number of true statements, check if it's consistent
    false_count = 7 - true_count
    
    # Initialize all statements as False
    statements = [False] * 7
    
    # Statement 1: At least 4 are true
    statements[0] = (true_count >= 4)
    
    # Statement 2: At most 1 is false
    statements[1] = (false_count <= 1)
    
    # Statement 3: Exactly 4 are true
    statements[2] = (true_count == 4)
    
    # Statement 4: Exactly 3 are false
    statements[3] = (false_count == 3)
    
    # Statement 5: Either Statement 3 or 4 is true, but not both
    statements[4] = (statements[2] != statements[3]) and (statements[2] or statements[3])
    
    # Statement 6: Number of true statements is prime
    statements[5] = is_prime(true_count)
    
    # Statement 7: Number of false statements is composite
    statements[6] = is_composite(false_count)
    
    # Check if the number of True statements matches true_count
    return sum(statements) == true_count

# Try all possible numbers of true statements (0 to 7)
solutions = []
for true_count in range(8):
    if check_consistency(true_count):
        solutions.append(true_count)
        print(f"Found solution with {true_count} true statements")

print(f"\nNumber of solutions: {len(solutions)}")
if solutions:
    print("Solutions found with these numbers of true statements:", solutions)