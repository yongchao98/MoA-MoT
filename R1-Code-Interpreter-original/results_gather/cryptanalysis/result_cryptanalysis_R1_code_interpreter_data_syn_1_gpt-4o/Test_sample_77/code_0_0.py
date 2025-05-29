from itertools import permutations

# Possible numbers and letters based on the conditions
possible_numbers = ['0', '1', '2']
possible_letters = [chr(i) for i in range(ord('D'), ord('Z') + 1) if i not in map(ord, 'WEKXFOACXY')]

# Function to check if a combination satisfies all conditions
def check_combination(numbers, letters):
    # Check each condition
    if numbers[0] in '65' or numbers[1] in '65':
        return False
    if numbers[0] not in '82' and numbers[1] not in '82':
        return False
    if numbers[0] not in '12' and numbers[1] not in '12':
        return False
    if numbers[0] in '97' or numbers[1] in '97':
        return False
    if numbers[0] in '38' or numbers[1] in '38':
        return False
    if numbers[0] not in '03' and numbers[1] not in '03':
        return False
    if numbers[0] in '85' or numbers[1] in '85':
        return False

    if letters[0] in 'WE' or letters[1] in 'WE':
        return False
    if letters[0] in 'RV' or letters[1] in 'RV':
        return False
    if letters[0] in 'KX' or letters[1] in 'KX':
        return False
    if letters[0] not in 'FO' and letters[1] not in 'FO':
        return False
    if letters[0] in 'AC' or letters[1] in 'AC':
        return False
    if letters[0] not in 'XY' and letters[1] not in 'XY':
        return False

    return True

# Try all permutations of the possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for let_perm in permutations(possible_letters, 2):
        if check_combination(num_perm, let_perm):
            print(list(num_perm + let_perm))
            break