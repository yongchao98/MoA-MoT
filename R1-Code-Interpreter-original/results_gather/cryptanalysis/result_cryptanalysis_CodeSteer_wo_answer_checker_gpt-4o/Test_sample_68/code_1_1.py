from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
conditions = [
    ('84KP', (0, 1, 1, 0)),  # (numbers correct, numbers wrong, letters correct, letters wrong)
    ('13GY', (1, 1, 0, 0)),
    ('31HM', (1, 1, 0, 0)),
    ('84RL', (0, 1, 1, 0)),
    ('02AN', (0, 0, 0, 0)),
    ('20BH', (0, 0, 0, 0)),
    ('17XY', (1, 1, 0, 0))
]

def check_combination(combo):
    for guess, feedback in conditions:
        num_correct = sum(c1 == c2 for c1, c2 in zip(combo[:2], guess[:2]))
        num_wrong = sum(c1 != c2 and c1 in guess[:2] for c1 in combo[:2])
        let_correct = sum(c1 == c2 for c1, c2 in zip(combo[2:], guess[2:]))
        let_wrong = sum(c1 != c2 and c1 in guess[2:] for c1 in combo[2:])
        
        if (num_correct, num_wrong, let_correct, let_wrong) != feedback:
            return False
    return True

# Generate all permutations of two numbers and two letters
for num_combo in permutations(numbers, 2):
    for let_combo in permutations(letters, 2):
        combo = num_combo + let_combo
        if check_combination(combo):
            print(list(combo))
            break