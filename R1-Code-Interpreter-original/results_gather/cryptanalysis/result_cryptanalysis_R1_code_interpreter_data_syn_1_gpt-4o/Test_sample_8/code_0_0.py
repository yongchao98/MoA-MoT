from itertools import permutations

# Possible numbers and letters based on the feedback
possible_numbers = [0, 1, 2]
possible_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'N']

# Function to check if a combination is valid
def is_valid_combination(numbers, letters):
    # Check against the first guess
    if numbers[0] >= 9 or numbers[1] >= 9:
        return False
    if letters[0] not in ['E', 'J'] and letters[1] not in ['E', 'J']:
        return False
    if letters[0] in ['E', 'J'] and letters[1] in ['E', 'J']:
        return False

    # Check against the second guess
    if numbers[0] >= 4 or numbers[1] >= 4:
        return False
    if letters[0] in ['O', 'M'] or letters[1] in ['O', 'M']:
        return False

    # Check against the third guess
    if numbers[0] not in [2, 0] and numbers[1] not in [2, 0]:
        return False
    if numbers[0] in [2, 0] and numbers[1] in [2, 0]:
        return False
    if letters[0] not in ['B', 'Q'] and letters[1] not in ['B', 'Q']:
        return False
    if letters[0] in ['B', 'Q'] and letters[1] in ['B', 'Q']:
        return False

    return True

# Iterate through all permutations of numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for let_perm in permutations(possible_letters, 2):
        if is_valid_combination(num_perm, let_perm):
            print(list(num_perm) + list(let_perm))
            break