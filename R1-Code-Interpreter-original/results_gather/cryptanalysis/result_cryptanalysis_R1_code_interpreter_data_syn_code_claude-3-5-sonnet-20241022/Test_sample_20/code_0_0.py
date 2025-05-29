def is_valid_guess(guess, target):
    # Check numbers
    nums_correct_pos = sum(1 for i in range(2) if guess[i] == target[i])
    nums_too_large = sum(1 for i in range(2) if int(guess[i]) > int(target[i]))
    nums_too_small = sum(1 for i in range(2) if int(guess[i]) < int(target[i]))
    
    # Check letters
    letters_correct_pos = sum(1 for i in range(2,4) if guess[i] == target[i])
    letters_too_late = sum(1 for i in range(2,4) if guess[i] > target[i])
    
    return nums_correct_pos, nums_too_large, nums_too_small, letters_correct_pos, letters_too_late

def check_conditions(candidate):
    # Condition 1: 35PQ
    c1 = is_valid_guess("35PQ", candidate)
    if not (c1[0] == 1 and c1[1] == 1 and c1[3] == 0 and c1[4] == 2):
        return False
    
    # Condition 2: 24XE
    c2 = is_valid_guess("24XE", candidate)
    if not (c2[0] == 1 and c2[2] == 1 and c2[3] == 1 and c2[4] == 1):
        return False
    
    # Condition 3: 96SQ
    c3 = is_valid_guess("96SQ", candidate)
    if not (c3[1] == 2 and c3[3] == 0 and c3[4] == 2):
        return False
    
    # Condition 4: 25XC
    c4 = is_valid_guess("25XC", candidate)
    if not (c4[0] == 0 and c4[3] == 1 and c4[4] == 1):
        return False
    
    return True

# Generate all possible combinations
import itertools

digits = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for nums in itertools.permutations(digits, 2):
    for lets in itertools.permutations(letters, 2):
        candidate = ''.join(nums + lets)
        if check_conditions(candidate):
            print(f"['{nums[0]}', '{nums[1]}', '{lets[0]}', '{lets[1]}']")