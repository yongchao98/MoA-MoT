def check_statements(true_count):
    # Returns whether each statement is true given the number of true statements
    false_count = 7 - true_count
    
    s1 = (true_count >= 7)  # At least 7 are true
    s2 = (false_count <= 6)  # At most 6 are false
    s3 = (true_count == 5)  # Exactly 5 are true
    s4 = (false_count == 6)  # Exactly 6 are false
    s5 = (s3 != s4)  # Either s3 or s4 is true, but not both
    s6 = true_count in [2, 3, 5, 7]  # Number of true statements is prime
    s7 = false_count in [4, 6, 8, 9]  # Number of false statements is composite
    
    statements = [s1, s2, s3, s4, s5, s6, s7]
    return statements

def count_solutions():
    solutions = []
    
    # Try all possible numbers of true statements (0 to 7)
    for true_count in range(8):
        statements = check_statements(true_count)
        actual_true_count = sum(1 for s in statements if s)
        
        # If the assumed number of true statements matches the actual count
        if actual_true_count == true_count:
            solutions.append(statements)
    
    print(f"Number of solutions: {len(solutions)}")
    if solutions:
        print("\nSolutions found:")
        for i, sol in enumerate(solutions, 1):
            print(f"Solution {i}: {sol}")
            print(f"Number of true statements: {sum(1 for s in sol if s)}")
            print(f"Number of false statements: {sum(1 for s in sol if not s)}")
            print()

count_solutions()