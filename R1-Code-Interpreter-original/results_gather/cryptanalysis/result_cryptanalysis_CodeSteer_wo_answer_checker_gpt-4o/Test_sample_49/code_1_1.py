from itertools import permutations

# Define possible numbers and letters
possible_numbers = ['0', '1', '4', '5']
possible_letters = ['U', 'A']

# Function to check if a combination satisfies the feedback
def satisfies_feedback(combo, guess, feedback):
    num_correct = 0
    num_wrong_position = 0
    num_too_large = 0
    num_too_small = 0
    letter_correct = 0
    letter_wrong = 0

    # Check numbers
    for i in range(2):
        if combo[i] == guess[i]:
            num_correct += 1
        elif combo[i] in guess[:2]:
            num_wrong_position += 1
        elif combo[i] > guess[i]:
            num_too_large += 1
        else:
            num_too_small += 1

    # Check letters
    for i in range(2, 4):
        if combo[i] == guess[i]:
            letter_correct += 1
        else:
            letter_wrong += 1

    # Compare with feedback
    if feedback == (num_correct, num_wrong_position, num_too_large, num_too_small, letter_correct, letter_wrong):
        return True
    return False

# Define feedback for each guess
feedbacks = [
    (0, 1, 1, 0, 0, 2),  # 75KB
    (0, 1, 0, 1, 2, 0),  # 15UA
    (0, 0, 0, 2, 0, 2),  # 09KJ
    (0, 1, 0, 1, 0, 2),  # 05PB
    (0, 0, 2, 0, 0, 2),  # 69DB
    (2, 0, 0, 0, 0, 2)   # 54ZF
]

# Iterate through all possible combinations
for num_combo in permutations(possible_numbers, 2):
    for letter_combo in permutations(possible_letters, 2):
        combo = num_combo + letter_combo
        if all(satisfies_feedback(combo, guess, feedback) for guess, feedback in zip(['75KB', '15UA', '09KJ', '05PB', '69DB', '54ZF'], feedbacks)):
            print(f"<<< {list(combo)} >>>")
            break