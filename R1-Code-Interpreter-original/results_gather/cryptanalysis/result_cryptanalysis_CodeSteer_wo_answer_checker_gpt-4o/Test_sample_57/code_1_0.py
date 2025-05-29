from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Conditions based on the feedback
conditions = [
    ('43HI', (1, 1, 0)),  # One number correct but wrong position, one number incorrect and too small, both letters incorrect
    ('45QH', (1, 1, 0)),  # One number correct but wrong position, one number incorrect, both letters incorrect
    ('01ZA', (0, 2, 1)),  # Both numbers incorrect and too small, one letter correct but wrong position, one letter incorrect and too late
    ('64RO', (1, 1, 0)),  # One number correct and correct position, one number incorrect, both letters incorrect and too late
    ('92SB', (1, 1, 0)),  # One number correct and correct position, one number incorrect and too small, both letters incorrect
    ('46DY', (1, 1, 0)),  # One number correct but wrong position, one number incorrect, both letters incorrect
    ('81CM', (0, 0, 1)),  # Both numbers incorrect, one letter correct and correct position, one letter incorrect
]

def check_condition(guess, condition):
    guess_numbers = guess[:2]
    guess_letters = guess[2:]
    correct_numbers = condition[0][:2]
    correct_letters = condition[0][2:]
    
    # Check numbers
    correct_number_position = sum(gn == cn for gn, cn in zip(guess_numbers, correct_numbers))
    correct_number_wrong_position = sum(gn in correct_numbers for gn in guess_numbers) - correct_number_position
    
    # Check letters
    correct_letter_position = sum(gl == cl for gl, cl in zip(guess_letters, correct_letters))
    correct_letter_wrong_position = sum(gl in correct_letters for gl in guess_letters) - correct_letter_position
    
    return (correct_number_position, correct_number_wrong_position, correct_letter_position) == condition[1]

# Iterate over all possible combinations of two numbers and two letters
for num_comb in permutations(numbers, 2):
    for letter_comb in permutations(letters, 2):
        candidate = ''.join(num_comb) + ''.join(letter_comb)
        if all(check_condition(candidate, condition) for condition in conditions):
            print(list(candidate))
            break