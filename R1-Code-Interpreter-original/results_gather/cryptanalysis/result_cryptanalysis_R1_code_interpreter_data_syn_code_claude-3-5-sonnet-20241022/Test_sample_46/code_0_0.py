def check_guess(actual, guess, feedback):
    # Unpack feedback
    num_correct_pos = sum(1 for i in range(2) if feedback[i] == 'correct_pos')
    num_correct_wrong_pos = sum(1 for i in range(2) if feedback[i] == 'correct_wrong_pos')
    num_incorrect = sum(1 for i in range(2) if feedback[i] == 'incorrect')
    
    letter_correct_wrong_pos = sum(1 for i in range(2,4) if feedback[i] == 'correct_wrong_pos')
    letter_too_early = sum(1 for i in range(2,4) if feedback[i] == 'too_early')
    letter_incorrect = sum(1 for i in range(2,4) if feedback[i] == 'incorrect')
    
    # Check numbers
    correct_pos = sum(1 for i in range(2) if actual[i] == guess[i])
    correct_wrong_pos = sum(1 for i in range(2) if guess[i] in actual[:2] and actual[i] != guess[i])
    incorrect = sum(1 for i in range(2) if guess[i] not in actual[:2])
    
    if (correct_pos != num_correct_pos or 
        correct_wrong_pos != num_correct_wrong_pos or 
        incorrect != num_incorrect):
        return False
    
    # Check letters
    l_correct_wrong_pos = sum(1 for i in range(2,4) if guess[i] in actual[2:] and actual[i] != guess[i])
    l_too_early = sum(1 for i in range(2,4) if ord(guess[i]) < ord(actual[i]))
    l_incorrect = sum(1 for i in range(2,4) if guess[i] not in actual[2:] and not ord(guess[i]) < ord(actual[i]))
    
    return (l_correct_wrong_pos == letter_correct_wrong_pos and 
            l_too_early == letter_too_early and 
            l_incorrect == letter_incorrect)

def is_valid_password(password):
    # Check guesses against known conditions
    guesses = [
        ('20BH', [('correct_wrong_pos', 'incorrect'), ('too_early', 'too_early')]),
        ('85IM', [('incorrect', 'incorrect'), ('too_early', 'too_early')]),
        ('09OQ', [('correct_pos', 'correct_pos'), ('correct_wrong_pos', 'incorrect')]),
        ('28BR', [('incorrect', 'incorrect'), ('correct_wrong_pos', 'too_early')])
    ]
    
    for guess, feedback in guesses:
        if not check_guess(password, guess, feedback):
            return False
    return True

# Generate all possible combinations
import string
valid_passwords = []
numbers = '0123456789'
letters = string.ascii_uppercase

for n1 in '09':  # We know these are correct from guess 3
    for n2 in '09':
        if n1 == n2:  # Numbers can't repeat
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:  # Letters can't repeat
                    continue
                password = [n1, n2, l1, l2]
                if is_valid_password(password):
                    valid_passwords.append(password)

print(valid_passwords)