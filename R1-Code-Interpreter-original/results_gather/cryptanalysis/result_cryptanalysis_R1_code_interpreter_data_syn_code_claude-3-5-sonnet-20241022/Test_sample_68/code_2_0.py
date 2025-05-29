from itertools import permutations
import string

def check_guess(guess_nums, guess_letters, actual_nums, actual_letters, num_feedback, letter_feedback):
    # Check numbers
    if num_feedback == "both incorrect":
        if any(n in actual_nums for n in guess_nums):
            return False
    elif num_feedback == "both incorrect and too small":
        if any(n in actual_nums for n in guess_nums) or any(int(n) >= min(int(x) for x in actual_nums) for n in guess_nums):
            return False
    elif num_feedback == "one correct in position":
        if not (guess_nums[0] == actual_nums[0] and guess_nums[1] != actual_nums[1]):
            return False
    elif num_feedback == "one correct wrong position":
        if not (guess_nums[1] == actual_nums[0] or guess_nums[0] == actual_nums[1]):
            return False
        if guess_nums[0] == actual_nums[0] or guess_nums[1] == actual_nums[1]:
            return False
    
    # Check letters
    if letter_feedback == "both incorrect":
        if any(l in actual_letters for l in guess_letters):
            return False
    elif letter_feedback == "one correct wrong position":
        correct_count = sum(l in actual_letters for l in guess_letters)
        if correct_count != 1:
            return False
        if guess_letters[0] == actual_letters[0] or guess_letters[1] == actual_letters[1]:
            return False
    
    return True

def is_valid_password(nums, letters):
    # Convert to lists for easier handling
    nums = list(nums)
    letters = list(letters)
    
    # Guess 1: 84KP
    if not check_guess(['8', '4'], ['K', 'P'], nums, letters, 
                      "both incorrect", "one correct wrong position"):
        return False
    
    # Guess 2: 13GY
    if not check_guess(['1', '3'], ['G', 'Y'], nums, letters,
                      "one correct in position", "both incorrect"):
        return False
    
    # Guess 3: 31HM
    if not check_guess(['3', '1'], ['H', 'M'], nums, letters,
                      "one correct wrong position", "both incorrect"):
        return False
    
    # Guess 4: 84RL
    if not check_guess(['8', '4'], ['R', 'L'], nums, letters,
                      "both incorrect", "one correct wrong position"):
        return False
    
    # Guess 5: 02AN
    if not (int(nums[0]) > 2 and int(nums[1]) > 2 and
            all(l not in letters for l in ['A', 'N'])):
        return False
    
    # Guess 6: 20BH
    if not (int(nums[0]) > 2 and int(nums[1]) > 2 and
            all(l not in letters for l in ['B', 'H'])):
        return False
    
    # Guess 7: 17XY
    if not check_guess(['1', '7'], ['X', 'Y'], nums, letters,
                      "one correct wrong position", "both incorrect"):
        return False
    
    # Additional constraints
    # P, L, X, Y should be too late in alphabet compared to actual letters
    if not all(ord(l) > max(ord(actual) for actual in letters) 
              for l in ['P', 'L', 'X', 'Y']):
        return False
    
    # B, H should be too early in alphabet
    if not all(ord(l) < min(ord(actual) for actual in letters) 
              for l in ['B', 'H']):
        return False
    
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = list(string.ascii_uppercase)

valid_passwords = []
for nums in permutations(numbers, 2):
    for lets in permutations(letters, 2):
        if is_valid_password(nums, lets):
            valid_passwords.append([nums[0], nums[1], lets[0], lets[1]])

print(valid_passwords)