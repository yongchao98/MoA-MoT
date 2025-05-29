from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

# Generate all possible combinations of two numbers and two letters
all_combinations = [
    (n1, n2, l1, l2)
    for n1, n2 in permutations(possible_numbers, 2)
    for l1, l2 in permutations(possible_letters, 2)
]

# Function to check constraints for each guess
def check_guess_1(n1, n2, l1, l2):
    return n1 not in '27' and n2 not in '27' and l1 not in 'JR' and l2 not in 'JR'

def check_guess_2(n1, n2, l1, l2):
    return n1 not in '67' and n2 not in '67' and ((l1 == 'G' and l2 != 'P') or (l2 == 'G' and l1 != 'P'))

def check_guess_3(n1, n2, l1, l2):
    return ((n1 == '0' and n2 != '9') or (n2 == '0' and n1 != '9')) and l1 not in 'CO' and l2 not in 'CO'

def check_guess_4(n1, n2, l1, l2):
    return n1 not in '57' and n2 not in '57' and l1 not in 'VJ' and l2 not in 'VJ'

def check_guess_5(n1, n2, l1, l2):
    return ((n1 == '0' and n2 != '6') or (n2 == '0' and n1 != '6')) and l1 not in 'HO' and l2 not in 'HO'

def check_guess_6(n1, n2, l1, l2):
    return n1 not in '67' and n2 not in '67' and ((l1 == 'Q' and l2 != 'J') or (l2 == 'Q' and l1 != 'J'))

def check_guess_7(n1, n2, l1, l2):
    return (n1 == '3' and n2 != '4') and l1 < 'T' and l2 < 'T'

# Filter combinations that satisfy all constraints
valid_combinations = [
    comb for comb in all_combinations
    if check_guess_1(*comb) and check_guess_2(*comb) and check_guess_3(*comb)
    and check_guess_4(*comb) and check_guess_5(*comb) and check_guess_6(*comb)
    and check_guess_7(*comb)
]

# Print the valid combination(s)
for combination in valid_combinations:
    print(f"<<< {list(combination)} >>>")