from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = ['0', '1', '2', '4']
possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1) if i not in map(ord, 'WEKXRVACFOY')]

# Feedback constraints
feedback = [
    ('65WE', (False, False, False, False)),
    ('82RV', (True, False, False, False)),
    ('12KX', (True, False, False, False)),
    ('97FO', (False, False, True, False)),
    ('38AC', (False, False, False, False)),
    ('03CX', (True, False, False, False)),
    ('85XY', (False, False, True, False))
]

def is_valid_combination(combination):
    for guess, (n1, n2, l1, l2) in feedback:
        # Check numbers
        if n1 and combination[0] != guess[0] and combination[0] != guess[1]:
            return False
        if n2 and combination[1] != guess[0] and combination[1] != guess[1]:
            return False
        # Check letters
        if l1 and combination[2] != guess[2] and combination[2] != guess[3]:
            return False
        if l2 and combination[3] != guess[2] and combination[3] != guess[3]:
            return False
    return True

# Generate all permutations of possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        combination = num_perm + letter_perm
        if is_valid_combination(combination):
            print(f"<<< {list(combination)} >>>")
            break