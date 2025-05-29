def check_guess(actual, guess, feedback):
    # Unpack feedback for numbers and letters separately
    num_feedback, letter_feedback = feedback
    
    # Count actual matches for numbers
    correct_pos = sum(1 for i in range(2) if actual[i] == guess[i])
    correct_wrong_pos = sum(1 for i in range(2) if guess[i] in actual[:2] and actual[i] != guess[i])
    incorrect = 2 - correct_pos - correct_wrong_pos
    
    # Compare number feedback
    if (correct_pos != num_feedback.count('correct_pos') or 
        correct_wrong_pos != num_feedback.count('correct_wrong_pos') or 
        incorrect != num_feedback.count('incorrect')):
        return False
    
    # Count actual matches for letters
    l_correct_wrong_pos = sum(1 for i in range(2,4) if guess[i] in ''.join(actual[2:]) and actual[i] != guess[i])
    l_too_early = sum(1 for i in range(2,4) if ord(guess[i]) < ord(actual[i]))
    l_incorrect = 2 - l_correct_wrong_pos - l_too_early
    
    # Compare letter feedback
    if (l_correct_wrong_pos != letter_feedback.count('correct_wrong_pos') or 
        l_too_early != letter_feedback.count('too_early') or 
        l_incorrect != letter_feedback.count('incorrect')):
        return False
    
    return True

def is_valid_combination(nums, letters):
    # Test against all guesses
    guesses = [
        ('20BH', (['correct_wrong_pos', 'incorrect'], ['too_early', 'too_early'])),
        ('85IM', (['incorrect', 'incorrect'], ['too_early', 'too_early'])),
        ('09OQ', (['correct_pos', 'correct_pos'], ['correct_wrong_pos', 'incorrect'])),
        ('28BR', (['incorrect', 'incorrect'], ['correct_wrong_pos', 'too_early']))
    ]
    
    password = list(nums + letters)
    for guess, feedback in guesses:
        if not check_guess(password, list(guess), feedback):
            return False
    return True

# We know from guess 3 that numbers are 09
numbers = '09'
valid_passwords = []

# Generate possible letter combinations
import string
letters = string.ascii_uppercase

# From the conditions:
# - Letters must be after B,H,I,M in alphabet
# - One letter must be O or Q (from guess 3)
# - One letter must be R (from guess 4)
# - Letters can't repeat

for l1 in letters:
    for l2 in letters:
        if l1 != l2 and ord(l1) > ord('M') and ord(l2) > ord('M'):
            if (l1 in 'OQ' or l2 in 'OQ') and (l1 == 'R' or l2 == 'R'):
                if is_valid_combination(numbers, l1 + l2):
                    valid_passwords.append(['0', '9', l1, l2])

print(valid_passwords)