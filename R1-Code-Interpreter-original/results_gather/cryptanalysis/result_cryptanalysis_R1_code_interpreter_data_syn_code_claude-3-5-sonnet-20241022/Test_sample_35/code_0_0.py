from itertools import permutations

def check_guess(guess, actual):
    # Check numbers
    numbers_correct_pos = sum(1 for i in range(2) if guess[i] == actual[i])
    numbers_correct = sum(1 for i in range(2) if str(guess[i]) in actual[:2]) - numbers_correct_pos
    
    # Check letters
    letters_correct_pos = sum(1 for i in range(2,4) if guess[i] == actual[i])
    letters_correct = sum(1 for i in range(2,4) if guess[i] in actual[2:]) - letters_correct_pos
    
    return numbers_correct_pos, numbers_correct, letters_correct_pos, letters_correct

def is_early_in_alphabet(letter, target):
    return ord(letter) < ord(target)

def validate_password(password):
    # Convert password to string format for comparison
    password_str = ''.join(str(x) for x in password)
    
    # Check guess 1: 14LZ
    if check_guess("14LZ", password_str) != (0, 2, 0, 0):
        return False
    
    # Check guess 2: 28EP
    pos_n, cor_n, pos_l, cor_l = check_guess("28EP", password_str)
    if not (pos_n == 0 and cor_n == 0 and pos_l == 1 and cor_l == 0):
        return False
    
    # Check guess 3: 14RW
    if check_guess("14RW", password_str) != (0, 2, 0, 0):
        return False
    
    # Check guess 4: 70XH
    pos_n, cor_n, pos_l, cor_l = check_guess("70XH", password_str)
    if not (pos_n == 0 and cor_n == 0 and pos_l == 1 and cor_l == 0):
        return False
    
    # Additional checks for letters being too early
    if 'E' in password_str[2:] or 'X' in password_str[2:]:
        return False
    
    return True

# Generate all possible combinations
numbers = [1, 4]  # We know these are the correct numbers
letters = [chr(i) for i in range(ord('A'), ord('Z')+1)]
letters = [l for l in letters if l not in 'LZRWEX']  # Remove known incorrect letters

valid_passwords = []
for num_perm in permutations(numbers, 2):
    for l1 in letters:
        for l2 in letters:
            if l1 != l2:  # Letters can't repeat
                password = [num_perm[0], num_perm[1], l1, l2]
                if validate_password(password):
                    valid_passwords.append(password)

print(valid_passwords)