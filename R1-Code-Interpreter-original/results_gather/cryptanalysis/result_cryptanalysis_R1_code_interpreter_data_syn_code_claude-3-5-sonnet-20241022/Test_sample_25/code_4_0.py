def super_detailed_check(solution):
    print(f"\nSuper detailed check for {solution}")
    nums = [int(solution[0]), int(solution[1])]
    letters = solution[2:]
    
    # Guess 1: 91SF
    print("\nGuess 1 (91SF):")
    print(f"- Neither 9 nor 1 in solution: {9 not in nums and 1 not in nums}")
    print(f"- S in correct position (pos 2): {solution[2] == 'S'}")
    print(f"- F not in solution: {'F' not in letters}")
    
    # Guess 2: 03HF
    print("\nGuess 2 (03HF):")
    if nums[0] == 0:
        print(f"- 0 correct in pos 0, 3 too small compared to {nums[1]}: {3 < nums[1]}")
    if nums[1] == 3:
        print(f"- 3 correct in pos 1, 0 too small compared to {nums[0]}: {0 < nums[0]}")
    print(f"- Neither H nor F in solution: {'H' not in letters and 'F' not in letters}")
    
    # Guess 3: 75CP
    print("\nGuess 3 (75CP):")
    print(f"- Neither 7 nor 5 in solution: {7 not in nums and 5 not in nums}")
    print(f"- C in solution but wrong position: {'C' in letters and solution[2] != 'C'}")
    print(f"- P not in solution: {'P' not in letters}")
    
    # Guess 4: 64YB
    print("\nGuess 4 (64YB):")
    print(f"- Either 6 or 4 correct in position: {nums[0] == 6 or nums[1] == 4}")
    print(f"- Other number not in solution at all")
    print(f"- Neither Y nor B in solution: {'Y' not in letters and 'B' not in letters}")

for sol in [['0', '4', 'S', 'C'], ['6', '3', 'S', 'C']]:
    super_detailed_check(sol)