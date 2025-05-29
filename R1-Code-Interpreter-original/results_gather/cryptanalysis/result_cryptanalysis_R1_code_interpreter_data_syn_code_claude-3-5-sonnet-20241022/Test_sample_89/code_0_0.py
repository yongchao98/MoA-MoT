from itertools import product

def check_guess(password, guess, feedback):
    # Unpack feedback
    num_correct = sum(1 for i in range(2) if password[i] == guess[i])
    let_correct = sum(1 for i in range(2,4) if password[i] == guess[i])
    let_wrong_pos = sum(1 for i in range(2,4) if guess[i] in password[2:] and guess[i] != password[i])
    
    if feedback == "both_num_incorrect":
        return num_correct == 0
    elif feedback == "one_num_correct_pos":
        return num_correct == 1
    elif feedback == "let_one_wrong_pos":
        return let_correct == 0 and let_wrong_pos == 1
    elif feedback == "both_let_incorrect":
        return let_correct == 0 and let_wrong_pos == 0
    return True

def is_valid_password(password):
    # Convert password to list of strings
    password = [str(x) for x in password]
    
    # Check all conditions
    if not check_guess(password, ['1','6','T','E'], "both_num_incorrect"): return False
    if not check_guess(password, ['7','2','Q','Y'], "one_num_correct_pos"): return False
    if not check_guess(password, ['2','3','A','X'], "one_num_correct_pos"): return False
    if not check_guess(password, ['7','1','K','Q'], "one_num_correct_pos"): return False
    if not check_guess(password, ['9','1','T','W'], "both_num_incorrect"): return False
    if not check_guess(password, ['5','1','P','L'], "both_num_incorrect"): return False
    
    # Additional specific conditions
    if password[0] == password[1]: return False  # numbers can't repeat
    if password[2] == password[3]: return False  # letters can't repeat
    
    # Check if one of Q,Y is in password (from guess 2)
    if not ('Q' in password[2:] or 'Y' in password[2:]): return False
    
    # Check if one of P,L is in password (from guess 6)
    if not ('P' in password[2:] or 'L' in password[2:]): return False
    
    return True

# Generate all possible combinations
numbers = range(0,10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
valid_passwords = []

for combo in product(numbers, numbers, letters, letters):
    if is_valid_password(combo):
        valid_passwords.append(combo)

print(valid_passwords)