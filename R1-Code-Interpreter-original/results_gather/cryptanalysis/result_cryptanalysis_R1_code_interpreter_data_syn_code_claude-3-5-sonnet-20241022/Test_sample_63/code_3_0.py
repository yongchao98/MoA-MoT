def test_guess(test_nums, test_letters, guess):
    nums, letters, conditions = guess
    
    # Count number conditions
    correct_pos_n = 0
    wrong_pos_n = 0
    too_large_n = 0
    too_small_n = 0
    
    # First digit
    if test_nums[0] == nums[0]:
        correct_pos_n += 1
    elif test_nums[0] == nums[1]:
        wrong_pos_n += 1
    elif test_nums[0] > nums[0]:
        too_large_n += 1
    else:
        too_small_n += 1
    
    # Second digit
    if test_nums[1] == nums[1]:
        correct_pos_n += 1
    elif test_nums[1] == nums[0]:
        wrong_pos_n += 1
    elif test_nums[1] > nums[1]:
        too_large_n += 1
    else:
        too_small_n += 1
    
    # Count letter conditions
    correct_pos_l = 0
    wrong_pos_l = 0
    too_early_l = 0
    too_late_l = 0
    
    # First letter
    if test_letters[0] == letters[0]:
        correct_pos_l += 1
    elif test_letters[0] == letters[1]:
        wrong_pos_l += 1
    elif ord(test_letters[0]) < ord(letters[0]):
        too_early_l += 1
    else:
        too_late_l += 1
    
    # Second letter
    if test_letters[1] == letters[1]:
        correct_pos_l += 1
    elif test_letters[1] == letters[0]:
        wrong_pos_l += 1
    elif ord(test_letters[1]) < ord(letters[1]):
        too_early_l += 1
    else:
        too_late_l += 1
    
    return (correct_pos_n == conditions[0][0] and 
            wrong_pos_n == conditions[0][1] and 
            too_large_n == conditions[0][2] and 
            too_small_n == conditions[0][3] and 
            correct_pos_l == conditions[1][0] and 
            wrong_pos_l == conditions[1][1] and 
            too_early_l == conditions[1][2] and 
            too_late_l == conditions[1][3])

# Define guesses with their conditions
guesses = [
    ((3,8), ('U','W'), ((0,1,1,0), (0,0,0,2))),  # 38UW
    ((6,9), ('O','G'), ((1,0,1,0), (1,0,1,0))),  # 69OG
    ((8,3), ('O','P'), ((1,0,1,0), (1,0,0,0))),  # 83OP
    ((4,8), ('S','P'), ((0,0,0,0), (0,0,0,0))),  # 48SP
    ((4,2), ('Z','E'), ((0,0,0,0), (0,0,0,0))),  # 42ZE
    ((0,1), ('W','G'), ((0,0,0,2), (0,0,0,0))),  # 01WG
    ((0,7), ('J','X'), ((0,0,0,0), (0,0,0,0))),  # 07JX
    ((3,1), ('H','K'), ((0,1,0,1), (0,0,2,0))),  # 31HK
    ((3,7), ('V','Y'), ((0,1,1,0), (0,0,0,2))),  # 37VY
    ((6,2), ('Q','H'), ((1,0,0,1), (0,0,0,0))),  # 62QH
    ((9,4), ('C','L'), ((0,0,0,0), (0,0,2,0))),  # 94CL
    ((7,1), ('T','E'), ((0,0,0,0), (0,1,1,0)))   # 71TE
]

# Key observations:
# 1. 6 must be in first position (69OG, 62QH both show one number correct in position)
# 2. 3 must be in second position (appears in wrong position in multiple guesses)
# 3. O must be in second position (69OG, 83OP both show one letter correct in position)
# 4. T must be in first position (appears in wrong position in 71TE)

# Test the deduced solution
solution = ([6,3], ['T','O'])
valid = True

print("Testing solution:", solution)
for guess in guesses:
    if not test_guess(solution[0], solution[1], guess):
        valid = False
        print(f"Failed at guess: {guess}")
        break

if valid:
    print("Solution verified!")
    print("Final answer:", [solution[0][0], solution[0][1], solution[1][0], solution[1][1]])