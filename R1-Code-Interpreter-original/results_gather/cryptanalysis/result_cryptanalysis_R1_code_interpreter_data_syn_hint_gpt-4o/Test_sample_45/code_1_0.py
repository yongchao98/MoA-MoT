from itertools import permutations

# Possible numbers and letters based on deductions
possible_numbers = ['0', '4']
possible_letters = ['S', 'U']

# Generate all permutations of the numbers and letters
number_permutations = permutations(possible_numbers)
letter_permutations = permutations(possible_letters)

# Define the conditions
conditions = [
    {"guess": "52OL", "feedback": (False, False, False, False)},
    {"guess": "20FP", "feedback": (True, False, False, False)},
    {"guess": "08EX", "feedback": (False, True, False, False)},
    {"guess": "58RN", "feedback": (False, False, False, False)},
    {"guess": "64UM", "feedback": (False, True, True, False)},
    {"guess": "02PV", "feedback": (False, True, False, False)},
    {"guess": "62CO", "feedback": (False, False, False, False)},
    {"guess": "91CD", "feedback": (False, False, False, False)},
    {"guess": "96SR", "feedback": (False, False, False, True)}
]

# Function to check if a guess satisfies the feedback
def check_guess(guess, feedback):
    return all(
        (guess[i] == correct_password[i]) == feedback[i]
        for i in range(4)
    )

# Iterate over all possible combinations
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        correct_password = list(num_perm) + list(let_perm)
        if all(check_guess(condition["guess"], condition["feedback"]) for condition in conditions):
            print(correct_password)
            break