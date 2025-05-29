def check_guess(guess, actual):
    feedback = {'numbers_correct_pos': 0, 'numbers_wrong_pos': 0, 'numbers_too_large': 0,
                'letters_correct_pos': 0, 'letters_wrong_pos': 0}
    
    # Check numbers (first two positions)
    for i in range(2):
        if guess[i] == actual[i]:
            feedback['numbers_correct_pos'] += 1
        elif guess[i] in actual[:2]:
            feedback['numbers_wrong_pos'] += 1
        elif int(guess[i]) > int(actual[i]):
            feedback['numbers_too_large'] += 1
    
    # Check letters (last two positions)
    for i in range(2, 4):
        if guess[i] == actual[i]:
            feedback['letters_correct_pos'] += 1
        elif guess[i] in actual[2:]:
            feedback['letters_wrong_pos'] += 1
    
    return feedback

def matches_conditions(combination):
    guesses = [
        ('07UC', {'numbers_wrong_pos': 1, 'numbers_too_large': 1, 'letters_correct_pos': 0, 'letters_wrong_pos': 0}),
        ('60FG', {'numbers_correct_pos': 1, 'numbers_too_large': 1, 'letters_correct_pos': 0, 'letters_wrong_pos': 0}),
        ('59VF', {'numbers_too_large': 2, 'letters_correct_pos': 1, 'letters_wrong_pos': 0}),
        ('85EK', {'numbers_too_large': 2, 'letters_correct_pos': 0, 'letters_wrong_pos': 1}),
        ('18MA', {'numbers_correct_pos': 1, 'numbers_too_large': 1, 'letters_correct_pos': 0, 'letters_wrong_pos': 0})
    ]
    
    for guess, expected in guesses:
        feedback = check_guess(guess, combination)
        # Compare relevant parts of feedback with expected
        if feedback['numbers_correct_pos'] != expected.get('numbers_correct_pos', 0) or \
           feedback['numbers_wrong_pos'] != expected.get('numbers_wrong_pos', 0) or \
           feedback['numbers_too_large'] != expected.get('numbers_too_large', 0) or \
           feedback['letters_correct_pos'] != expected.get('letters_correct_pos', 0) or \
           feedback['letters_wrong_pos'] != expected.get('letters_wrong_pos', 0):
            return False
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
valid_combinations = []

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # Numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # Letters can't repeat
                        combination = n1 + n2 + l1 + l2
                        if matches_conditions(combination):
                            valid_combinations.append([n1, n2, l1, l2])

print(valid_combinations)