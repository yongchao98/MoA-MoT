def verify_solution(combination):
    # Count true statements
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Helper function to check if a number is prime
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0:
                return False
        return True
    
    # Helper function to check if a number is composite
    def is_composite(n):
        return n > 1 and not is_prime(n)
    
    # Evaluate each statement
    s1 = true_count >= 6
    s2 = false_count <= 1
    s3 = true_count == 6
    s4 = false_count == 5
    s5 = (combination[2] != combination[3])
    s6 = is_prime(true_count)
    s7 = is_composite(false_count)
    
    evaluations = [s1, s2, s3, s4, s5, s6, s7]
    
    print(f"\nFor combination {combination}:")
    print(f"True count: {true_count}, False count: {false_count}")
    print("Statement evaluations:")
    print(f"1. At least 6 true: {s1}")
    print(f"2. At most 1 false: {s2}")
    print(f"3. Exactly 6 true: {s3}")
    print(f"4. Exactly 5 false: {s4}")
    print(f"5. Statement 3 XOR 4: {s5}")
    print(f"6. True count is prime: {s6}")
    print(f"7. False count is composite: {s7}")
    
    return combination == [1 if s else 0 for s in evaluations]

# Check both solutions
solutions = [[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1]]
valid_solutions = [s for s in solutions if verify_solution(s)]
print(f"\nNumber of valid solutions: {len(valid_solutions)}")