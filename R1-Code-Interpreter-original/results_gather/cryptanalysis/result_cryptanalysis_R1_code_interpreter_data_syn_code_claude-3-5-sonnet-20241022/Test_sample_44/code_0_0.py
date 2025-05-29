import itertools

def check_numbers_against_guess(guess_nums, actual_nums):
    correct_pos = sum(1 for i in range(2) if guess_nums[i] == actual_nums[i])
    correct_wrong_pos = sum(1 for x in guess_nums if x in actual_nums) - correct_pos
    return correct_pos, correct_wrong_pos

def check_letters_against_guess(guess_letters, actual_letters):
    correct_pos = sum(1 for i in range(2) if guess_letters[i] == actual_letters[i])
    correct_wrong_pos = sum(1 for x in guess_letters if x in actual_letters) - correct_pos
    return correct_pos, correct_wrong_pos

def is_valid_password(password):
    # Convert password to numbers and letters
    nums = password[:2]
    letters = password[2:]
    
    # Check all conditions
    # Condition 1: 71HT
    if check_numbers_against_guess(['7','1'], nums) != (1,0): return False
    if int(nums[0]) <= 1 or int(nums[1]) <= 1: return False
    
    # Condition 2: 19ZS
    if any(x in nums for x in ['1','9']): return False
    
    # Condition 6: 76KR (most informative)
    if nums != ['7','6']: return False
    if check_letters_against_guess(['K','R'], letters) != (0,1): return False
    
    # Condition 12: 65FG
    if check_letters_against_guess(['F','G'], letters) != (1,0): return False
    
    # Condition 11: 26SR
    if check_letters_against_guess(['S','R'], letters) != (0,1): return False
    
    return True

# Generate all possible combinations
numbers = ['7','6']  # We know these from condition 6
letters = [chr(i) for i in range(65, 91)]  # A-Z
valid_passwords = []

for l1, l2 in itertools.permutations(letters, 2):
    password = numbers + [l1, l2]
    if is_valid_password(password):
        valid_passwords.append(password)

print(valid_passwords)