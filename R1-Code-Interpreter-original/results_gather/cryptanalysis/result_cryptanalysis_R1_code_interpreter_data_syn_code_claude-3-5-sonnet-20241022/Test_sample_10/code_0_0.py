from itertools import permutations

def check_guess(password, guess, feedback):
    # Unpack feedback
    num_correct = feedback['numbers']
    num_position = feedback.get('num_position', [])
    num_range = feedback.get('num_range', [])
    let_correct = feedback['letters']
    let_position = feedback.get('let_position', [])
    let_range = feedback.get('let_range', [])
    
    # Check numbers
    correct_nums = 0
    for i in range(2):
        if guess[i] == password[i]:
            correct_nums += 1
        elif guess[i] in password[:2]:
            if 'correct_wrong_pos' in num_position:
                correct_nums += 1
        elif 'too_small' in num_range and int(guess[i]) < int(password[i]):
            continue
        elif 'too_large' in num_range and int(guess[i]) > int(password[i]):
            continue
    if correct_nums != num_correct:
        return False
    
    # Check letters
    correct_lets = 0
    for i in range(2,4):
        if guess[i] == password[i]:
            correct_lets += 1
        elif guess[i] in password[2:]:
            if 'correct_wrong_pos' in let_position:
                correct_lets += 1
        elif 'too_late' in let_range and guess[i] > password[i]:
            continue
    if correct_lets != let_correct:
        return False
    
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Store guesses and their feedback
guesses = [
    ('25KS', {'numbers': 0, 'letters': 0}),
    ('26BP', {'numbers': 0, 'letters': 1, 'let_position': ['correct_wrong_pos'], 'let_range': ['too_late']}),
    ('39ON', {'numbers': 1, 'letters': 1, 'num_position': ['correct_wrong_pos'], 'num_range': ['too_small'], 
              'let_position': ['correct_wrong_pos'], 'let_range': ['too_late']})
]

possible_passwords = []
for nums in permutations(numbers, 2):
    for lets in permutations(letters, 2):
        password = nums + lets
        valid = True
        for guess, feedback in guesses:
            if not check_guess(password, guess, feedback):
                valid = False
                break
        if valid:
            possible_passwords.append(password)

print(possible_passwords)