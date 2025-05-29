from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = ['0', '1', '2', '3']
possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'N']

# Function to check if a combination satisfies all conditions
def is_valid_combination(numbers, letters):
    # Check against first guess: 97EJ
    if numbers[0] >= '9' or numbers[1] >= '7':
        return False
    if letters[0] == 'E' or letters[1] == 'J':
        return False
    if 'E' not in letters and 'J' not in letters:
        return False

    # Check against second guess: 43OM
    if numbers[0] >= '4' or numbers[1] >= '3':
        return False
    if letters[0] >= 'O' or letters[1] >= 'M':
        return False

    # Check against third guess: 20BQ
    if numbers[0] == '2' and numbers[1] == '0':
        return False
    if 'B' in letters or 'Q' in letters:
        return False
    if 'E' not in letters and 'J' not in letters:
        return False

    return True

# Iterate through all permutations of two numbers and two letters
for num_perm in permutations(possible_numbers, 2):
    for let_perm in permutations(possible_letters, 2):
        if is_valid_combination(num_perm, let_perm):
            password = list(num_perm + let_perm)
            print(f"<<< {password} >>>")
            break