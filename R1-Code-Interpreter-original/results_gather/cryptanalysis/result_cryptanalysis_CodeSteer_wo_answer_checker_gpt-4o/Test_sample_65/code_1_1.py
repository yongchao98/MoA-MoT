from itertools import permutations

# Possible numbers and letters based on the feedback
possible_numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

# Feedback conditions
conditions = [
    ("95OK", "both numbers are correct but in the wrong positions; both letters are incorrect and too early in the alphabet."),
    ("70EZ", "both numbers are incorrect; both letters are incorrect."),
    ("95RN", "both numbers are correct but in the wrong positions; both letters are incorrect and too early in the alphabet."),
    ("08EF", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("16OY", "both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet."),
    ("18GK", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("19SZ", "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect."),
    ("63EK", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("43CR", "both numbers are incorrect and too small; both letters are incorrect and too early in the alphabet."),
    ("41FV", "both numbers are incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet.")
]

# Function to check if a combination satisfies all conditions
def satisfies_conditions(combo):
    num1, num2, let1, let2 = combo
    for guess, feedback in conditions:
        g_num1, g_num2, g_let1, g_let2 = guess
        if "both numbers are correct but in the wrong positions" in feedback:
            if not (num1 == g_num2 and num2 == g_num1):
                return False
        if "both numbers are incorrect" in feedback:
            if num1 in guess or num2 in guess:
                return False
        if "one number is correct and in the correct position" in feedback:
            if not (num1 == g_num1 or num2 == g_num2):
                return False
        if "one number is incorrect and too small" in feedback:
            if not (num1 > g_num1 or num2 > g_num2):
                return False
        if "both letters are incorrect and too early in the alphabet" in feedback:
            if let1 in guess or let2 in guess:
                return False
        if "one letter is correct and in the correct position" in feedback:
            if not (let1 == g_let1 or let2 == g_let2):
                return False
        if "one letter is correct but in the wrong position" in feedback:
            if not (let1 == g_let2 or let2 == g_let1):
                return False
    return True

# Generate all possible combinations of two numbers and two letters
for num_combo in permutations(possible_numbers, 2):
    for let_combo in permutations(possible_letters, 2):
        combo = num_combo + let_combo
        if satisfies_conditions(combo):
            print(f"<<< {list(combo)} >>>")
            break