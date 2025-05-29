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

# Define the correct combination
correct_combination = "31CN"

# Function to check if a guess satisfies the feedback
def check_guess(guess, feedback, correct_combination):
    correct_numbers = sum(1 for i in range(2) if guess[i] == correct_combination[i])
    correct_letters = sum(1 for i in range(2, 4) if guess[i] == correct_combination[i])
    return (correct_numbers, correct_letters) == feedback

# Verify all conditions
all_conditions_satisfied = all(check_guess(cond["guess"], cond["feedback"], correct_combination) for cond in conditions)

print(all_conditions_satisfied)