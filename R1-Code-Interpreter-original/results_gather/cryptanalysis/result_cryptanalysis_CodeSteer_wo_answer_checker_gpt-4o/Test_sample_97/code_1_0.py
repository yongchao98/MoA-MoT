from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback from each guess
feedback = [
    ('71BQ', (1, 1, 0, 0)),  # (correct number in wrong position, incorrect number too small, incorrect letters too early)
    ('29BQ', (1, 1, 0, 0)),
    ('69BE', (1, 1, 0, 0)),
    ('40QT', (0, 2, 1, 0)),
    ('48NP', (0, 0, 0, 0)),
    ('15FQ', (0, 2, 0, 0)),
    ('83PF', (0, 0, 0, 0)),
    ('96MJ', (1, 1, 0, 0)),
    ('32YJ', (0, 2, 0, 0)),
    ('16SU', (0, 2, 0, 0)),
    ('89XJ', (1, 1, 0, 0)),
    ('36XJ', (0, 2, 0, 0)),
    ('35ZU', (0, 2, 0, 0)),
    ('36FV', (0, 2, 0, 0))
]

# Function to check if a combination is valid
def is_valid_combination(combo, feedback):
    for guess, (num_correct_wrong_pos, num_incorrect_small, letter_correct_wrong_pos, letter_incorrect_early) in feedback:
        num_correct = sum(c1 == c2 for c1, c2 in zip(combo[:2], guess[:2]))
        num_wrong_pos = sum(c in guess[:2] for c in combo[:2]) - num_correct
        num_small = sum(c < g for c, g in zip(combo[:2], guess[:2]))
        
        letter_correct = sum(c1 == c2 for c1, c2 in zip(combo[2:], guess[2:]))
        letter_wrong_pos = sum(c in guess[2:] for c in combo[2:]) - letter_correct
        letter_early = sum(c < g for c, g in zip(combo[2:], guess[2:]))
        
        if not (num_wrong_pos == num_correct_wrong_pos and num_small == num_incorrect_small and
                letter_wrong_pos == letter_correct_wrong_pos and letter_early == letter_incorrect_early):
            return False
    return True

# Iterate over all possible combinations
for num_combo in permutations(numbers, 2):
    for letter_combo in permutations(letters, 2):
        combo = num_combo + letter_combo
        if is_valid_combination(combo, feedback):
            print(f"<<< {list(combo)} >>>")
            break