def detailed_verify(solution):
    print(f"\nDetailed verification for {solution}")
    
    # Guess 1: 91SF
    # - both numbers incorrect
    # - one letter correct and in correct position
    # - one letter incorrect
    nums = [int(solution[0]), int(solution[1])]
    print("\nGuess 1 (91SF):")
    print(f"Numbers 9,1 not in solution: {9 not in nums and 1 not in nums}")
    print(f"S in correct position: {'S' == solution[2]}")
    print(f"F not in solution: {'F' not in solution[2:]}")

    # Guess 2: 03HF
    # - one number correct and in correct position
    # - one number incorrect and too small
    # - both letters incorrect
    print("\nGuess 2 (03HF):")
    print(f"Either 0 or 3 in correct position: {0 == nums[0] or 3 == nums[1]}")
    print(f"One number too small compared to solution")
    print(f"H and F not in solution: {'H' not in solution[2:] and 'F' not in solution[2:]}")

    # Guess 3: 75CP
    # - both numbers incorrect
    # - one letter correct but wrong position
    # - one letter incorrect
    print("\nGuess 3 (75CP):")
    print(f"Numbers 7,5 not in solution: {7 not in nums and 5 not in nums}")
    print(f"Either C or P in solution but wrong position: {('C' in solution[2:] and solution.index('C') != 2) or ('P' in solution[2:] and solution.index('P') != 3)}")

    # Guess 4: 64YB
    # - one number correct and in correct position
    # - one number incorrect
    # - both letters incorrect
    print("\nGuess 4 (64YB):")
    print(f"Either 6 or 4 in correct position: {6 == nums[0] or 4 == nums[1]}")
    print(f"Y and B not in solution: {'Y' not in solution[2:] and 'B' not in solution[2:]}")

    # Additional check for guess 2
    if 0 in nums:
        print("\nCritical check for guess 2:")
        print(f"If 0 is in solution, the other number in 03 must be too small")
    if 3 in nums:
        print(f"If 3 is in solution, the other number in 03 must be too small")

print("Testing first solution:")
detailed_verify(['0', '4', 'S', 'C'])
print("\nTesting second solution:")
detailed_verify(['6', '3', 'S', 'C'])