def final_check(solution):
    nums = [int(solution[0]), int(solution[1])]
    print(f"\nFinal verification for {solution}")
    
    # Guess 2 (03HF) detailed analysis
    print("\nDetailed analysis of guess 2 (03HF):")
    # Check if exactly one number is correct and in position
    correct_pos = (nums[0] == 0) + (nums[1] == 3)
    print(f"Numbers in solution: {nums}")
    print(f"Number of correct positions (should be 1): {correct_pos}")
    
    # Check the "too small" condition
    if nums[0] == 0:  # If 0 is correct
        print(f"3 compared to {nums[1]}: {'too small' if 3 < nums[1] else 'NOT too small'}")
    elif nums[1] == 3:  # If 3 is correct
        print(f"0 compared to {nums[0]}: {'too small' if 0 < nums[0] else 'NOT too small'}")
    else:
        print("Neither 0 nor 3 is in correct position!")

    # Guess 4 (64YB) detailed analysis
    print("\nDetailed analysis of guess 4 (64YB):")
    correct_pos = (nums[0] == 6) + (nums[1] == 4)
    print(f"Number of correct positions (should be 1): {correct_pos}")
    print(f"6 in solution: {6 in nums}")
    print(f"4 in solution: {4 in nums}")

solutions = [['0', '4', 'S', 'C'], ['6', '3', 'S', 'C']]
for sol in solutions:
    final_check(sol)