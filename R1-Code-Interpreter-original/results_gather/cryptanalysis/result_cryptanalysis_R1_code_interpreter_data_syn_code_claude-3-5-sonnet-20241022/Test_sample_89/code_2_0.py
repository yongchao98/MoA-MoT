def check_letter_position(guess_letters, password_letters, correct_wrong_pos):
    # Count letters in wrong position
    wrong_pos_count = sum(1 for i, letter in enumerate(guess_letters) 
                         for j, pass_letter in enumerate(password_letters) 
                         if i != j and letter == pass_letter)
    return wrong_pos_count == correct_wrong_pos

def is_valid_combination(password):
    nums, letters = password[:2], password[2:]
    
    # Check guess 2 (72QY)
    if not (nums[0] == '7' and  # one number correct in position
            ('Q' in letters or 'Y' in letters) and  # one letter in wrong position
            check_letter_position(['Q', 'Y'], letters, 1)):  # exactly one letter in wrong pos
        return False
    
    # Check guess 6 (51PL)
    if not (('P' in letters or 'L' in letters) and  # one letter must be P or L
            check_letter_position(['P', 'L'], letters, 1)):  # exactly one letter in wrong pos
        return False
    
    # Check guess 3 (23AX)
    if nums[1] != '3':  # second number must be 3
        return False
    
    # Check guess 4 (71KQ)
    if nums[0] != '7':  # first number must be 7
        return False
    
    # Additional check for letter that's too early in alphabet from guess 6
    if 'P' in letters:
        other_letter = letters[0] if letters[1] == 'P' else letters[1]
        if other_letter < 'P':
            return False
    if 'L' in letters:
        other_letter = letters[0] if letters[1] == 'L' else letters[1]
        if other_letter < 'L':
            return False
            
    return True

# Test all possible combinations
numbers = ['7', '3']  # We know these are fixed
letters = ['L', 'P', 'Q', 'Y']
valid_passwords = []

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # letters can't repeat
            password = ['7', '3', l1, l2]
            if is_valid_combination(password):
                valid_passwords.append(password)

print(valid_passwords)