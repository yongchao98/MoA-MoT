def verify_combination(statements):
    # statements is a list of 7 booleans
    true_count = sum(statements)
    false_count = 7 - true_count
    
    # Check each statement's claim against the actual combination
    valid = True
    
    # Statement 1: At least 4 are true
    valid = valid and (statements[0] == (true_count >= 4))
    
    # Statement 2: At most 0 are false
    valid = valid and (statements[1] == (false_count <= 0))
    
    # Statement 3: Exactly 4 are true
    valid = valid and (statements[2] == (true_count == 4))
    
    # Statement 4: Exactly 3 are false
    valid = valid and (statements[3] == (false_count == 3))
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    valid = valid and (statements[4] == ((statements[2] or statements[3]) and not (statements[2] and statements[3])))
    
    # Statement 6: Number of true statements is prime
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    valid = valid and (statements[5] == is_prime(true_count))
    
    # Statement 7: Number of false statements is composite
    def is_composite(n):
        return n > 1 and not is_prime(n)
    valid = valid and (statements[6] == is_composite(false_count))
    
    return valid

# Try all possible combinations
solutions = []
from itertools import product
for combo in product([False, True], repeat=7):
    if verify_combination(combo):
        solutions.append((sum(combo), combo))

print(f"Number of solutions: {len(solutions)}")
for true_count, combo in solutions:
    print(f"\nSolution with {true_count} true statements:")
    for i, value in enumerate(combo, 1):
        print(f"Statement {i}: {'True' if value else 'False'}")