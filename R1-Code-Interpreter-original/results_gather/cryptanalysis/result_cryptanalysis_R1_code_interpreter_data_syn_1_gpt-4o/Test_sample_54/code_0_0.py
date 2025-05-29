from itertools import permutations

# Possible numbers and letters based on the analysis
possible_numbers = [1, 2, 3, 8, 9]
possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']

# Function to check if a guess matches the feedback
def check_guess(guess, feedback):
    numbers, letters = guess[:2], guess[2:]
    num_feedback, let_feedback = feedback

    # Check numbers
    if num_feedback == "both incorrect":
        if any(n in numbers for n in possible_numbers):
            return False
    elif num_feedback == "one correct":
        if not any(n in numbers for n in possible_numbers):
            return False

    # Check letters
    if let_feedback == "both incorrect":
        if any(l in letters for l in possible_letters):
            return False
    elif let_feedback == "one correct":
        if not any(l in letters for l in possible_letters):
            return False

    return True

# Feedback for each guess
feedbacks = [
    ("both incorrect", "both incorrect"),
    ("both incorrect", "one correct"),
    ("both incorrect", "both incorrect"),
    ("both incorrect", "one correct"),
    ("one correct", "one correct"),
    ("both incorrect", "both incorrect"),
    ("one correct", "both incorrect")
]

# All possible combinations of numbers and letters
all_combinations = permutations(possible_numbers, 2)
all_letter_combinations = permutations(possible_letters, 2)

# Find the correct combination
for num_comb in all_combinations:
    for let_comb in all_letter_combinations:
        guess = num_comb + let_comb
        if all(check_guess(guess, feedback) for feedback in feedbacks):
            print(guess)
            break