from itertools import permutations

# Possible numbers and letters
possible_numbers = set(range(10))
possible_letters = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

# Constraints based on feedback
# 1. Guess: 85UJ
possible_numbers -= {8, 5}
possible_letters -= {'U', 'J'}

# 2. Guess: 19BP
possible_numbers -= {1, 9}
possible_letters -= {'P'}
correct_letter_wrong_position = {'B'}

# 3. Guess: 78LU
possible_numbers -= {8}
correct_number_wrong_position = {7}
correct_letter_correct_position = {'L'}
possible_letters -= {'U'}

# 4. Guess: 61VZ
possible_numbers -= {6}
possible_letters -= {'V', 'Z'}

# 5. Guess: 41EL
possible_numbers -= {4}
possible_letters -= {'E'}
correct_letter_wrong_position.add('L')

# 6. Guess: 12PR
possible_numbers -= {2}
possible_letters -= {'P', 'R'}

# Remaining possible numbers and letters
remaining_numbers = list(possible_numbers)
remaining_letters = list(possible_letters)

# Function to check if a combination satisfies all constraints
def is_valid_combination(numbers, letters):
    # Check numbers
    if 7 in numbers and numbers.index(7) == 0:
        return False
    # Check letters
    if 'B' in letters and letters.index('B') == 1:
        return False
    if 'L' not in letters or letters.index('L') != 1:
        return False
    return True

# Generate all permutations of the remaining numbers and letters
for num_perm in permutations(remaining_numbers, 2):
    for let_perm in permutations(remaining_letters, 2):
        if is_valid_combination(num_perm, let_perm):
            password = list(map(str, num_perm)) + list(let_perm)
            print(f"<<< {password} >>>")
            break