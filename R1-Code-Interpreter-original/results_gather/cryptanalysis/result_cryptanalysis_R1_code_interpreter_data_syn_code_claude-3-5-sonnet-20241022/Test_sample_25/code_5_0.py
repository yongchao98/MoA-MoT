def final_crucial_check(solution):
    print(f"\nCrucial check for {solution}")
    nums = [int(solution[0]), int(solution[1])]
    
    # For Guess 4 (64YB):
    # - one number correct and in correct position
    # - one number incorrect and incorrect (meaning it shouldn't appear anywhere)
    
    print("\nGuess 4 (64YB) analysis:")
    print(f"Numbers in solution: {nums}")
    
    # Check 6
    six_correct_pos = (nums[0] == 6)
    six_anywhere = 6 in nums
    
    # Check 4
    four_correct_pos = (nums[1] == 4)
    four_anywhere = 4 in nums
    
    print(f"6 in correct position: {six_correct_pos}")
    print(f"6 appears anywhere: {six_anywhere}")
    print(f"4 in correct position: {four_correct_pos}")
    print(f"4 appears anywhere: {four_anywhere}")
    
    # The crucial test:
    # If one number is correct in position, the other number should not appear anywhere
    if six_correct_pos:
        print(f"If 6 is correct in position, 4 should not appear anywhere: {not four_anywhere}")
    if four_correct_pos:
        print(f"If 4 is correct in position, 6 should not appear anywhere: {not six_anywhere}")

solutions = [['0', '4', 'S', 'C'], ['6', '3', 'S', 'C']]
for sol in solutions:
    final_crucial_check(sol)