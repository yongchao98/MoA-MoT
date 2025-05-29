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

# Define the correct password
correct_password = ["4", "0", "S", "U"]

# Function to check if a guess satisfies the feedback
def check_guess(guess, feedback):
    return all(
        (guess[i] == correct_password[i]) == feedback[i]
        for i in range(4)
    )

# Check all conditions
all_conditions_satisfied = all(
    check_guess(condition["guess"], condition["feedback"])
    for condition in conditions
)

print(all_conditions_satisfied)