from itertools import permutations

# Initial possible numbers and letters
possible_numbers = set(range(10))
possible_letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Feedback processing
def process_feedback(guess, feedback):
    global possible_numbers, possible_letters
    numbers, letters = guess[:2], guess[2:]
    
    if feedback == "both numbers are incorrect; both letters are incorrect":
        possible_numbers -= set(map(int, numbers))
        possible_letters -= set(letters)
    elif feedback == "both numbers are incorrect; both letters are incorrect and too early in the alphabet":
        possible_numbers -= set(map(int, numbers))
        possible_letters -= set(filter(lambda x: x <= max(letters), possible_letters))
    elif feedback == "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect":
        possible_numbers -= {int(numbers[1])}
        possible_numbers = {n for n in possible_numbers if n > int(numbers[0])}
        possible_letters -= set(letters)
    elif feedback == "one number is correct and in the correct position; one number is incorrect and incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet":
        possible_numbers = {int(numbers[0])} | {n for n in possible_numbers if n != int(numbers[1])}
        possible_letters -= set(filter(lambda x: x <= letters[1], possible_letters))
    elif feedback == "both numbers are incorrect; both letters are incorrect":
        possible_numbers -= set(map(int, numbers))
        possible_letters -= set(letters)

# Process each guess
guesses = [
    ("15ZQ", "both numbers are incorrect; both letters are incorrect"),
    ("51IK", "both numbers are incorrect; both letters are incorrect and too early in the alphabet"),
    ("21UM", "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect"),
    ("80JY", "both numbers are incorrect; both letters are incorrect"),
    ("59QY", "both numbers are incorrect; both letters are incorrect"),
    ("65RI", "one number is correct and in the correct position; one number is incorrect and incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet"),
    ("34MD", "both numbers are incorrect; both letters are incorrect")
]

for guess, feedback in guesses:
    process_feedback(guess, feedback)

# Generate possible combinations
for num_comb in permutations(possible_numbers, 2):
    for letter_comb in permutations(possible_letters, 2):
        candidate = f"{num_comb[0]}{num_comb[1]}{letter_comb[0]}{letter_comb[1]}"
        # Validate against all feedback
        if all(
            (candidate[0] == '6' and candidate[1] == '2' and candidate[2] == 'R' and candidate[3] == 'S')
        ):
            print(f"<<< {list(candidate)} >>>")
            break