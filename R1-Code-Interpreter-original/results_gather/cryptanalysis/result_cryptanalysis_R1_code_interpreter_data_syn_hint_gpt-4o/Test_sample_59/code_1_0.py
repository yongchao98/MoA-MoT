from itertools import permutations

# Define the conditions
conditions = [
    {"guess": "20IX", "feedback": (0, 0)},
    {"guess": "09YT", "feedback": (1, 0)},
    {"guess": "78LK", "feedback": (0, 0)},
    {"guess": "15LD", "feedback": (1, 0)},
    {"guess": "03CO", "feedback": (0, 1)},
    {"guess": "21XS", "feedback": (1, 0)},
    {"guess": "48TP", "feedback": (0, 0)},
    {"guess": "25LH", "feedback": (0, 0)},
    {"guess": "72NO", "feedback": (0, 1)},
    {"guess": "65RX", "feedback": (0, 0)}
]

# Possible numbers and letters based on deductions
possible_numbers = ['1', '3']
possible_letters = ['C', 'N']

# Function to check if a guess satisfies the feedback
def check_guess(guess, feedback, correct_combination):
    correct_numbers = sum(1 for i in range(2) if guess[i] == correct_combination[i])
    correct_letters = sum(1 for i in range(2, 4) if guess[i] == correct_combination[i])
    return (correct_numbers, correct_letters) == feedback

# Generate all permutations of the possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        combination = ''.join(num_perm) + ''.join(letter_perm)
        if all(check_guess(cond["guess"], cond["feedback"], combination) for cond in conditions):
            print(combination)