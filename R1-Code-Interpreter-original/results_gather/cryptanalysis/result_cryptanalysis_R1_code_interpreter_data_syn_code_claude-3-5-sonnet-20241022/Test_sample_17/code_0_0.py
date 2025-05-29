def check_guess(guess, actual):
    feedback = {
        'correct_num_wrong_pos': 0,
        'incorrect_num_too_small': 0,
        'incorrect_num_too_large': 0,
        'correct_letter_right_pos': 0,
        'incorrect_letter_too_early': 0,
        'incorrect_letter': 0
    }
    
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            feedback['correct_num_wrong_pos'] += 1  # This is wrong, should be correct position
        elif int(guess[i]) in [int(actual[0]), int(actual[1])]:
            feedback['correct_num_wrong_pos'] += 1
        elif int(guess[i]) < min(int(actual[0]), int(actual[1])):
            feedback['incorrect_num_too_small'] += 1
        else:
            feedback['incorrect_num_too_large'] += 1
            
    # Check letters
    for i in range(2, 4):
        if guess[i] == actual[i]:
            feedback['correct_letter_right_pos'] += 1
        elif guess[i] < actual[i]:
            feedback['incorrect_letter_too_early'] += 1
        else:
            feedback['incorrect_letter'] += 1
            
    return feedback

def matches_conditions(password):
    # Convert password list to proper format for checking
    password_str = ''.join(str(x) for x in password)
    
    # Check all conditions
    guess1 = check_guess("10PL", password_str)
    if not (guess1['correct_num_wrong_pos'] == 1 and 
            guess1['incorrect_num_too_small'] == 1 and 
            guess1['incorrect_letter_too_early'] == 2):
        return False
        
    guess2 = check_guess("98QX", password_str)
    if not (guess2['incorrect_num_too_large'] == 2 and 
            guess2['correct_letter_right_pos'] == 1 and 
            guess2['incorrect_letter'] == 1):
        return False
        
    guess3 = check_guess("10KL", password_str)
    if not (guess3['correct_num_wrong_pos'] == 1 and 
            guess3['incorrect_num_too_small'] == 1 and 
            guess3['incorrect_letter_too_early'] == 2):
        return False
        
    guess4 = check_guess("16NY", password_str)
    if not (guess4['correct_num_wrong_pos'] == 2 and 
            guess4['incorrect_letter'] == 2):
        return False
        
    return True

# Generate all possible combinations
import itertools
import string

numbers = ['1', '6']
letters = list(string.ascii_uppercase)  # All uppercase letters

# Try all possible combinations
for num_perm in itertools.permutations(numbers, 2):
    for letter_combo in itertools.combinations(letters, 2):
        for letter_perm in itertools.permutations(letter_combo, 2):
            password = list(num_perm + letter_perm)
            if matches_conditions(password):
                print(password)