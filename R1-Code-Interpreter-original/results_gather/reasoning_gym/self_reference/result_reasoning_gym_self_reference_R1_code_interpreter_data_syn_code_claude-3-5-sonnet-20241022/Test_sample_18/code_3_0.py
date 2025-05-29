def verify_both_solutions():
    # Solution 1: All false
    print("Verifying Solution 1 (all false):")
    true_count_1 = 0
    false_count_1 = 7
    print(f"Statement 1 should be false: {not (true_count_1 >= 4)}")
    print(f"Statement 2 should be false: {not (false_count_1 <= 1)}")
    print(f"Statement 3 should be false: {not (true_count_1 == 4)}")
    print(f"Statement 4 should be false: {not (false_count_1 == 3)}")
    print(f"Statement 5 should be false: {not (False)}")  # since neither 3 nor 4 is true
    print(f"Statement 6 should be false: {not (true_count_1 in {2,3,5,7})}")
    print(f"Statement 7 should be false: {not (false_count_1 in {4,6})}")
    
    print("\nVerifying Solution 2 (only statement 7 true):")
    true_count_2 = 1
    false_count_2 = 6
    print(f"Statement 1 should be false: {not (true_count_2 >= 4)}")
    print(f"Statement 2 should be false: {not (false_count_2 <= 1)}")
    print(f"Statement 3 should be false: {not (true_count_2 == 4)}")
    print(f"Statement 4 should be false: {not (false_count_2 == 3)}")
    print(f"Statement 5 should be false: {not (False)}")  # since neither 3 nor 4 is true
    print(f"Statement 6 should be false: {not (true_count_2 in {2,3,5,7})}")
    print(f"Statement 7 should be true: {false_count_2 in {4,6}}")

verify_both_solutions()