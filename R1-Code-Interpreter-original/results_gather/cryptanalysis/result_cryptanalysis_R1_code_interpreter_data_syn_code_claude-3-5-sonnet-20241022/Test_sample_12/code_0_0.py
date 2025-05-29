from itertools import product
import string

def check_guess(candidate, guess, number_feedback, letter_feedback):
    # Extract numbers and letters from both
    cand_nums = candidate[:2]
    cand_lets = candidate[2:]
    guess_nums = guess[:2]
    guess_lets = guess[2:]
    
    # Check numbers
    correct_nums = 0
    correct_pos_nums = 0
    for i in range(2):
        if guess_nums[i] in cand_nums:
            correct_nums += 1
            if guess_nums[i] == cand_nums[i]:
                correct_pos_nums += 1
    
    # Check letters
    correct_lets = 0
    correct_pos_lets = 0
    for i in range(2):
        if guess_lets[i] in cand_lets:
            correct_lets += 1
            if guess_lets[i] == cand_lets[i]:
                correct_pos_lets += 1
    
    return (correct_nums, correct_pos_nums, correct_lets, correct_pos_lets)

def is_valid_candidate(candidate):
    # Check all conditions
    
    # Guess 1: 86WU
    nums, pos_nums, lets, pos_lets = check_guess(candidate, "86WU")
    if nums != 0 or lets != 1 or pos_lets != 0:
        return False
    if not ('W' in candidate[2:] or 'U' in candidate[2:]):
        return False
        
    # Guess 2: 94MF
    nums, pos_nums, lets, pos_lets = check_guess(candidate, "94MF")
    if nums != 1 or pos_nums != 0 or lets != 0:
        return False
    if ord(candidate[2]) <= ord('M') or ord(candidate[3]) <= ord('F'):
        return False
        
    # Guess 3: 80UV
    nums, pos_nums, lets, pos_lets = check_guess(candidate, "80UV")
    if nums != 0 or lets != 0:
        return False
        
    # Guess 4: 47TH
    nums, pos_nums, lets, pos_lets = check_guess(candidate, "47TH")
    if nums != 0 or lets != 1 or pos_lets != 1:
        return False
        
    # Guess 5: 03SU
    nums, pos_nums, lets, pos_lets = check_guess(candidate, "03SU")
    if nums != 1 or pos_nums != 0 or lets != 0:
        return False
    if int(candidate[0]) <= 0 or int(candidate[1]) <= 3:
        return False
    
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = string.ascii_uppercase

valid_candidates = []
for n1, n2 in product(numbers, numbers):
    if n1 != n2:  # No repeating numbers
        for l1, l2 in product(letters, letters):
            if l1 != l2:  # No repeating letters
                candidate = n1 + n2 + l1 + l2
                if is_valid_candidate(candidate):
                    valid_candidates.append([n1, n2, l1, l2])

print(valid_candidates)