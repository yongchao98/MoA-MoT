from itertools import permutations

# Define possible numbers and letters
possible_numbers = '0123456789'
possible_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
conditions = [
    ("71HT", (1, 1, 0, 0)),
    ("19ZS", (0, 0, 0, 0)),
    ("48PA", (0, 0, 0, 0)),
    ("49FO", (0, 0, 0, 0)),
    ("25IP", (0, 0, 0, 0)),
    ("76KR", (2, 0, 1, 0)),
    ("17BY", (0, 1, 0, 0)),
    ("31ZP", (0, 0, 0, 0)),
    ("02NE", (0, 0, 0, 0)),
    ("36ZF", (1, 1, 0, 0)),
    ("26SR", (1, 1, 1, 0)),
    ("65FG", (0, 1, 1, 0))
]

# Function to check if a guess satisfies the feedback
def satisfies_conditions(guess, feedback):
    num_correct_pos = sum(guess[i] == feedback[i] for i in range(4))
    num_correct = sum(guess[i] in feedback for i in range(4)) - num_correct_pos
    return (num_correct_pos, num_correct)

# Iterate through all possible combinations
for num1, num2 in permutations(possible_numbers, 2):
    for let1, let2 in permutations(possible_letters, 2):
        guess = num1 + num2 + let1 + let2
        if all(satisfies_conditions(guess, cond[0]) == cond[1] for cond in conditions):
            print(f"<<< [{num1}, {num2}, {let1}, {let2}] >>>")
            break