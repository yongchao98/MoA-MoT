from itertools import permutations

# Define the constraints based on the feedback
def satisfies_constraints(guess, feedback, candidate):
    for i, (g, f, c) in enumerate(zip(guess, feedback, candidate)):
        if f == 'correct':
            if g != c:
                return False
        elif f == 'wrong position':
            if g == c or g not in candidate:
                return False
        elif f == 'incorrect':
            if g in candidate:
                return False
    return True

# Define the feedback for each guess
feedbacks = [
    ('79AG', ['wrong position', 'incorrect', 'incorrect', 'incorrect']),
    ('15BN', ['correct', 'incorrect', 'incorrect', 'incorrect']),
    ('46QN', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('16BC', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('81BA', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('69IH', ['wrong position', 'incorrect', 'incorrect', 'incorrect']),
    ('05WO', ['correct', 'incorrect', 'incorrect', 'incorrect']),
    ('74PV', ['incorrect', 'incorrect', 'correct', 'incorrect']),
    ('26BF', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('63TH', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('74HT', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('63OY', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('06XR', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('63HS', ['incorrect', 'incorrect', 'incorrect', 'incorrect']),
    ('13EX', ['incorrect', 'incorrect', 'correct', 'incorrect']),
]

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Generate all possible combinations of two numbers and two letters
for num_comb in permutations(numbers, 2):
    for let_comb in permutations(letters, 2):
        candidate = num_comb + let_comb
        if all(satisfies_constraints(guess, feedback, candidate) for guess, feedback in feedbacks):
            print(f"<<< {list(candidate)} >>>")
            break