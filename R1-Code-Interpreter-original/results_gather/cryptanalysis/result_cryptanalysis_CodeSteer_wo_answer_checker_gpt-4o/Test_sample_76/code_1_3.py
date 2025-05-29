from itertools import permutations

# Define the constraints based on the feedback
def satisfies_constraints(guess, feedback, candidate):
    correct_number_position = feedback[0]
    correct_number_wrong_position = feedback[1]
    correct_letter_position = feedback[2]
    correct_letter_wrong_position = feedback[3]

    # Split guess into numbers and letters
    guess_numbers = guess[:2]
    guess_letters = guess[2:]

    # Split candidate into numbers and letters
    candidate_numbers = candidate[:2]
    candidate_letters = candidate[2:]

    # Check numbers
    numbers_correct_position = sum(1 for i in range(2) if candidate_numbers[i] == guess_numbers[i])
    numbers_wrong_position = sum(1 for i in range(2) if candidate_numbers[i] in guess_numbers and candidate_numbers[i] != guess_numbers[i])

    # Check letters
    letters_correct_position = sum(1 for i in range(2) if candidate_letters[i] == guess_letters[i])
    letters_wrong_position = sum(1 for i in range(2) if candidate_letters[i] in guess_letters and candidate_letters[i] != guess_letters[i])

    # Adjust for double-counting in wrong position checks
    numbers_wrong_position -= numbers_correct_position
    letters_wrong_position -= letters_correct_position

    return (numbers_correct_position == correct_number_position and
            numbers_wrong_position == correct_number_wrong_position and
            letters_correct_position == correct_letter_position and
            letters_wrong_position == correct_letter_wrong_position)

# Define the feedback for each guess
guesses = [
    ("05XU", (0, 0, 0, 0)),
    ("80UC", (0, 1, 0, 0)),
    ("61HP", (0, 1, 0, 0)),
    ("38JO", (1, 0, 0, 0)),
    ("41TO", (0, 1, 0, 0)),
    ("39PF", (0, 0, 0, 1)),
    ("51KR", (0, 1, 0, 0)),
    ("13BE", (1, 0, 0, 0)),
    ("83AZ", (0, 1, 0, 0)),
    ("54LV", (0, 0, 1, 0))
]

# Possible numbers and letters based on constraints
possible_numbers = set(range(10)) - {0, 5, 6, 3, 4, 9}
possible_letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ") - set("XUCPHTOJBEAZKR")

# Try all permutations of the possible numbers and letters
for numbers in permutations(possible_numbers, 2):
    for letters in permutations(possible_letters, 2):
        candidate = list(numbers) + list(letters)
        if all(satisfies_constraints(guess, feedback, candidate) for guess, feedback in guesses):
            print(f"<<< {candidate} >>>")
            break