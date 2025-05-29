from itertools import permutations

# Possible numbers and letters based on the conditions
possible_numbers = [0, 4, 5, 9]
possible_letters = ['E', 'F', 'H', 'I', 'M', 'O', 'Q', 'R', 'S', 'U']

# All possible permutations of numbers and letters
number_permutations = list(permutations(possible_numbers, 2))
letter_permutations = list(permutations(possible_letters, 2))

# Check the conditions
def check_combination(numbers, letters):
    conditions = [
        (87, 'VX', False, False),
        (26, 'AP', False, False),
        (98, 'ZT', True, False),
        (98, 'HL', True, True),
        (93, 'BN', True, False),
        (37, 'CG', False, False),
        (17, 'KW', False, False),
        (92, 'XC', True, False),
        (30, 'NX', False, False),
        (74, 'DJ', True, False)
    ]
    
    for num, let, num_correct, let_correct in conditions:
        num_str = str(num)
        if num_correct:
            if not ((numbers[0] == int(num_str[0]) or numbers[1] == int(num_str[1]))):
                return False
        else:
            if numbers[0] == int(num_str[0]) or numbers[1] == int(num_str[1]):
                return False
        
        if let_correct:
            if not ((letters[0] == let[0] or letters[1] == let[1])):
                return False
        else:
            if letters[0] == let[0] or letters[1] == let[1]:
                return False
    
    return True

# Find the correct combination
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        if check_combination(num_perm, let_perm):
            print(num_perm, let_perm)
            break