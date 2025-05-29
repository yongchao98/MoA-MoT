def check_statements(true_count):
    false_count = 7 - true_count
    
    # Statement 1: At least 4 are true
    s1 = true_count >= 4
    
    # Statement 2: At most 0 are false
    s2 = false_count <= 0
    
    # Statement 3: Exactly 4 are true
    s3 = true_count == 4
    
    # Statement 4: Exactly 3 are false
    s4 = false_count == 3
    
    # Statement 5: Either Statement 3 or Statement 4 is true, but not both
    s5 = (s3 or s4) and not (s3 and s4)
    
    # Statement 6: Number of true statements is prime
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    s6 = is_prime(true_count)
    
    # Statement 7: Number of false statements is composite
    def is_composite(n):
        return n > 1 and not is_prime(n)
    s7 = is_composite(false_count)
    
    # Count how many statements are actually true
    actual_true = sum([s1, s2, s3, s4, s5, s6, s7])
    
    # Check if this combination is consistent
    return actual_true == true_count

# Try all possible numbers of true statements (0 to 7)
solutions = []
for true_count in range(8):
    if check_statements(true_count):
        solutions.append(true_count)

print(f"Number of solutions: {len(solutions)}")
print(f"Solutions (number of true statements): {solutions}")