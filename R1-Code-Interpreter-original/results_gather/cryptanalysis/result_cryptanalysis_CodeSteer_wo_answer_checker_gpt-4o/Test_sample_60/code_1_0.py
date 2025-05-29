from itertools import permutations

# Define the constraints based on the feedback
def satisfies_constraints(numbers, letters):
    # Constraint 1: 15PU
    if not ((numbers[0] == '1' or numbers[1] == '1') and numbers[0] != '5' and numbers[1] != '5' and letters[0] < 'P' and letters[1] < 'U'):
        return False
    # Constraint 2: 43EH
    if not (numbers[0] < '4' and numbers[1] < '3' and ('E' in letters or 'H' in letters) and (letters[0] > 'E' or letters[1] > 'E')):
        return False
    # Constraint 3: 94HG
    if not (numbers[0] < '9' and numbers[1] < '4' and letters[1] == 'H' and (letters[0] > 'G' or letters[1] > 'G')):
        return False
    # Constraint 4: 36NW
    if not (numbers[0] < '3' and numbers[1] < '6' and letters[0] < 'N' and letters[1] < 'W'):
        return False
    # Constraint 5: 39FU
    if not (numbers[0] < '3' and numbers[1] < '9' and 'F' not in letters and 'U' not in letters):
        return False
    # Constraint 6: 27JZ
    if not (numbers[0] < '2' and numbers[1] < '7' and ('J' in letters) and (letters[0] < 'Z' and letters[1] < 'Z')):
        return False
    return True

# Possible numbers and letters
possible_numbers = ['0', '1']
possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Generate all permutations of two numbers and two letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        if satisfies_constraints(num_perm, letter_perm):
            password = list(num_perm) + list(letter_perm)
            print(password)
            break