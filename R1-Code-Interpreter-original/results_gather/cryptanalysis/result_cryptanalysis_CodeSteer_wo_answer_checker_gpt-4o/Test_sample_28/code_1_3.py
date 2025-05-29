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

# Validation function
def validate_candidate(candidate, guesses):
    for guess, feedback in guesses:
        numbers, letters = guess[:2], guess[2:]
        candidate_numbers, candidate_letters = candidate[:2], candidate[2:]
        
        if feedback == "both numbers are incorrect; both letters are incorrect":
            if any(n in candidate_numbers for n in numbers) or any(l in candidate_letters for l in letters):
                return False
        elif feedback == "both numbers are incorrect; both letters are incorrect and too early in the alphabet":
            if any(n in candidate_numbers for n in numbers) or any(l in candidate_letters for l in letters):
                return False
            if any(l <= max(letters) for l in candidate_letters):
                return False
        elif feedback == "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect":
            if candidate_numbers[1] == numbers[0] or candidate_numbers[0] <= numbers[0]:
                return False
            if any(l in candidate_letters for l in letters):
                return False
        elif feedback == "one number is correct and in the correct position; one number is incorrect and incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet":
            if candidate_numbers[0] != numbers[0] or candidate_numbers[1] == numbers[1]:
                return False
            if candidate_letters[1] == letters[0] or candidate_letters[0] <= letters[1]:
                return False
        elif feedback == "both numbers are incorrect; both letters are incorrect":
            if any(n in candidate_numbers for n in numbers) or any(l in candidate_letters for l in letters):
                return False
    return True

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
        if validate_candidate(candidate, guesses):
            print(f"<<< {list(candidate)} >>>")
            break