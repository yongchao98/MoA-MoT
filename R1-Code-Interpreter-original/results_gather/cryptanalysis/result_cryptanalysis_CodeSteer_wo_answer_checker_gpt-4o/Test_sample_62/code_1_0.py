from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = [0, 6, 8]
possible_letters = ['O', 'I']

# Constraints based on feedback
def is_valid_combination(numbers, letters):
    # Check each guess against the feedback
    # Guess 1: 13GO
    if numbers[0] in [1, 3] or numbers[1] in [1, 3] or letters[0] != 'O' or letters[1] == 'G':
        return False
    # Guess 2: 26LC
    if numbers[0] == 2 or numbers[1] == 6 or 'L' in letters or 'C' in letters:
        return False
    # Guess 3: 95SX
    if numbers[0] in [9, 5] or numbers[1] in [9, 5] or 'S' in letters or 'X' in letters:
        return False
    # Guess 4: 94UG
    if numbers[0] in [9, 4] or numbers[1] in [9, 4] or 'U' in letters or 'G' in letters:
        return False
    # Guess 5: 12AT
    if numbers[0] in [1, 2] or numbers[1] in [1, 2] or 'A' in letters or 'T' in letters:
        return False
    # Guess 6: 61KJ
    if numbers[0] != 6 or numbers[1] == 1 or 'K' in letters or 'J' in letters:
        return False
    # Guess 7: 56GL
    if numbers[0] == 5 or numbers[1] == 6 or 'G' in letters or 'L' in letters:
        return False
    # Guess 8: 37BM
    if numbers[0] in [3, 7] or numbers[1] in [3, 7] or 'B' in letters or 'M' in letters:
        return False
    # Guess 9: 24IW
    if numbers[0] in [2, 4] or numbers[1] in [2, 4] or letters[0] != 'I' or letters[1] == 'W':
        return False
    return True

# Generate all permutations of the possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        if is_valid_combination(num_perm, letter_perm):
            password = list(num_perm) + list(letter_perm)
            print(f"<<< {password} >>>")
            break