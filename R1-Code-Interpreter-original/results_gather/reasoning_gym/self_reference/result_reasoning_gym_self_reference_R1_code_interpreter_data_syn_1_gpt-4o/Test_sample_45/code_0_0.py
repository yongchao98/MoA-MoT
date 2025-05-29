# Check the validity of each configuration
def check_statements(true_count):
    # Calculate false count
    false_count = 7 - true_count
    
    # Evaluate each statement
    s1 = true_count >= 5
    s2 = false_count <= 1
    s3 = true_count == 3
    s4 = false_count == 1
    s5 = (s3 or s4) and not (s3 and s4)
    s6 = true_count in [2, 3, 5, 7]  # Prime numbers
    s7 = false_count in [4, 6, 8]  # Composite numbers
    
    # Count true statements
    true_statements = sum([s1, s2, s3, s4, s5, s6, s7])
    
    # Check if the number of true statements matches the true_count
    return true_statements == true_count

# Check all possible true counts
possible_solutions = sum(check_statements(tc) for tc in range(8))
print(possible_solutions)