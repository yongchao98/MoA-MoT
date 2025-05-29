from itertools import permutations

# Possible numbers and letters based on feedback
possible_numbers = [0, 1, 3, 4, 7, 8]  # Excluding 2, 5, 6, 9
possible_letters = [chr(i) for i in range(65, 91) if i not in (75, 83, 80)]  # Excluding K, S, P

# Function to check if a combination is valid based on feedback
def is_valid_combination(numbers, letters):
    # Check against the first guess: 25KS
    if set(numbers).intersection({2, 5}) or set(letters).intersection({'K', 'S'}):
        return False
    
    # Check against the second guess: 26BP
    if set(numbers).intersection({2, 6}):
        return False
    if 'B' in letters and letters.index('B') == 1:
        return False
    if 'P' in letters:
        return False
    
    # Check against the third guess: 39ON
    if 3 in numbers and numbers.index(3) == 0:
        return False
    if 9 in numbers:
        return False
    if 'O' in letters and letters.index('O') == 2:
        return False
    if 'N' in letters:
        return False
    
    return True

# Iterate through all permutations of two numbers and two letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        if is_valid_combination(num_perm, letter_perm):
            password = list(num_perm) + list(letter_perm)
            print(f"<<< {password} >>>")
            break