from itertools import permutations

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
        ("38KE", (0, 0, 0, 0)),  # both numbers wrong, both letters wrong and early
        ("23BW", (0, 0, 0, 0)),  # both numbers wrong, both letters wrong
        ("60YQ", (0, 1, 0, 0)),  # one number wrong pos, one wrong, both letters wrong
        ("38LE", (0, 0, 0, 0)),  # both numbers wrong, both letters wrong and early
        ("18KS", (0, 0, 0, 0)),  # both numbers wrong, both letters wrong
        ("40NX", (0, 1, 1, 0)),  # one number wrong pos, one wrong, one letter correct pos, one wrong
        ("21JH", (0, 0, 0, 0)),  # both numbers wrong, both letters wrong and early
        ("73UR", (0, 1, 0, 1))   # one number wrong pos, one wrong, one letter wrong pos, one wrong
    ]
    
    for guess, feedback in guesses:
        result = check_guess(guess, combination)
        if result != feedback:
            return False
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

valid_combinations = []

# Try all possible combinations
for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        combination = n1 + n2 + l1 + l2
                        if is_consistent(combination):
                            valid_combinations.append([n1, n2, l1, l2])

print(valid_combinations)