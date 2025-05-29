from itertools import permutations

def is_letter_too_early(letter):
    early_threshold = 'L'  # Based on feedback from guesses
    return letter <= early_threshold

def is_letter_too_late(letter):
    late_threshold = 'W'  # Based on feedback
    return letter > late_threshold

def check_guess(guess, actual):
    # Returns feedback similar to the given conditions
    numbers_correct = 0
    numbers_wrong_pos = 0
    letters_correct = 0
    letters_wrong_pos = 0
    
    # Check numbers (first two positions)
    for i in range(2):
        if guess[i] == actual[i]:
            numbers_correct += 1
        elif guess[i] in actual[:2]:
            numbers_wrong_pos += 1
            
    # Check letters (last two positions)
    for i in range(2, 4):
        if guess[i] == actual[i]:
            letters_correct += 1
        elif guess[i] in actual[2:]:
            letters_wrong_pos += 1
            
    return (numbers_correct, numbers_wrong_pos, letters_correct, letters_wrong_pos)

def is_consistent(combination):
    # Check against all known guesses
    guesses = [
        ("38KE", (0, 0, 0, 0)),
        ("23BW", (0, 0, 0, 0)),
        ("60YQ", (0, 1, 0, 0)),
        ("38LE", (0, 0, 0, 0)),
        ("18KS", (0, 0, 0, 0)),
        ("40NX", (0, 1, 1, 0)),
        ("21JH", (0, 0, 0, 0)),
        ("73UR", (0, 1, 0, 1))
    ]
    
    for guess, feedback in guesses:
        result = check_guess(guess, combination)
        if result != feedback:
            return False
            
        # Additional letter position constraints
        if guess == "40NX":
            if combination[3] > 'W':  # X was too late
                return False
        if guess in ["38KE", "38LE", "21JH"]:
            if not (combination[2] > 'L' and combination[3] > 'L'):  # Letters should be later than L
                return False
    
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'MNOPQRSTUVW'  # Restricted letter range based on constraints

valid_combinations = []

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        combination = n1 + n2 + l1 + l2
                        if is_consistent(combination):
                            valid_combinations.append([n1, n2, l1, l2])

print(valid_combinations)