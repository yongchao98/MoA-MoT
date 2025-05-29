def test_combination(nums, letters, guess_data):
    guess_nums = guess_data[0]
    guess_letters = guess_data[1]
    feedback_nums = guess_data[2][0]
    feedback_letters = guess_data[2][1]
    
    # Test numbers
    correct_pos_n = 0
    wrong_pos_n = 0
    too_large_n = 0
    too_small_n = 0
    
    # Check first number
    if nums[0] == guess_nums[0]:
        correct_pos_n += 1
    elif nums[0] == guess_nums[1]:
        wrong_pos_n += 1
    elif nums[0] > guess_nums[0]:
        too_large_n += 1
    else:
        too_small_n += 1
        
    # Check second number
    if nums[1] == guess_nums[1]:
        correct_pos_n += 1
    elif nums[1] == guess_nums[0]:
        wrong_pos_n += 1
    elif nums[1] > guess_nums[1]:
        too_large_n += 1
    else:
        too_small_n += 1
    
    # Test letters
    correct_pos_l = 0
    wrong_pos_l = 0
    too_early_l = 0
    too_late_l = 0
    
    # Check first letter
    if letters[0] == guess_letters[0]:
        correct_pos_l += 1
    elif letters[0] == guess_letters[1]:
        wrong_pos_l += 1
    elif ord(letters[0]) < ord(guess_letters[0]):
        too_early_l += 1
    else:
        too_late_l += 1
        
    # Check second letter
    if letters[1] == guess_letters[1]:
        correct_pos_l += 1
    elif letters[1] == guess_letters[0]:
        wrong_pos_l += 1
    elif ord(letters[1]) < ord(guess_letters[1]):
        too_early_l += 1
    else:
        too_late_l += 1
    
    # Compare with feedback
    return (correct_pos_n == feedback_nums[0] and 
            wrong_pos_n == feedback_nums[1] and 
            too_large_n == feedback_nums[2] and 
            too_small_n == feedback_nums[3] and 
            correct_pos_l == feedback_letters[0] and 
            wrong_pos_l == feedback_letters[1] and 
            too_early_l == feedback_letters[2] and 
            too_late_l == feedback_letters[3])

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

# Test all possible combinations that match our key observations
possible_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2:
                    continue
                
                valid = True
                for guess in guesses:
                    if not test_combination([n1,n2], [l1,l2], guess):
                        valid = False
                        break
                
                if valid:
                    possible_solutions.append([n1, n2, l1, l2])
                    print(f"Found valid solution: {[n1, n2, l1, l2]}")

if not possible_solutions:
    print("No solutions found")
else:
    print(f"Total solutions found: {len(possible_solutions)}")