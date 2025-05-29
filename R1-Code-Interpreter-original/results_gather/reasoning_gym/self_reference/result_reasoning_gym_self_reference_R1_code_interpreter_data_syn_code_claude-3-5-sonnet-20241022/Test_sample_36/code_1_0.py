def evaluate_and_explain(true_count):
    false_count = 7 - true_count
    
    # Evaluate each statement
    s1 = true_count >= 4
    s2 = false_count <= 7
    s3 = true_count == 4
    s4 = false_count == 3
    s5 = (s3 or s4) and not (s3 and s4)
    
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n ** 0.5) + 1):
            if n % i == 0: return False
        return True
    
    s6 = is_prime(true_count)
    s7 = not is_prime(false_count) and false_count > 1
    
    statements = [s1, s2, s3, s4, s5, s6, s7]
    actual_true = sum(statements)
    
    if actual_true == true_count:
        print(f"\nFor {true_count} true statements:")
        print(f"Statement 1 (≥4 true): {s1}")
        print(f"Statement 2 (≤7 false): {s2}")
        print(f"Statement 3 (=4 true): {s3}")
        print(f"Statement 4 (=3 false): {s4}")
        print(f"Statement 5 (3 XOR 4): {s5}")
        print(f"Statement 6 (true count prime): {s6}")
        print(f"Statement 7 (false count composite): {s7}")
        print(f"Total true statements: {actual_true}")
        return True
    return False

for true_count in range(8):
    evaluate_and_explain(true_count)