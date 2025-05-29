def analyze_solution(true_count):
    false_count = 7 - true_count
    
    # For each statement, determine if it's true or false
    s1 = (true_count >= 4)  # "At least 4 are true"
    s2 = (false_count <= 1)  # "At most 1 is false"
    s3 = (true_count == 4)  # "Exactly 4 are true"
    s4 = (false_count == 3)  # "Exactly 3 are false"
    s5 = (s3 != s4 and (s3 or s4))  # "Either Statement 3 or 4 is true, but not both"
    s6 = true_count in {2, 3, 5, 7}  # "Number of true statements is prime"
    s7 = false_count in {4, 6, 8, 9}  # "Number of false statements is composite"
    
    statements = [s1, s2, s3, s4, s5, s6, s7]
    actual_true_count = sum(statements)
    
    print(f"\nAnalyzing solution with {true_count} assumed true statements:")
    print(f"Statement values: {statements}")
    print(f"Actual number of true statements: {actual_true_count}")
    print(f"Is consistent: {actual_true_count == true_count}")
    
    return actual_true_count == true_count

# Analyze both solutions
for true_count in [0, 1]:
    analyze_solution(true_count)