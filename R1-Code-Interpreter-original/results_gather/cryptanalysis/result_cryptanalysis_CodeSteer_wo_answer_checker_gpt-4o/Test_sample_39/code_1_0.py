from itertools import permutations

def is_valid_combination(numbers, letters):
    # Check against the first guess: 87ID
    # One number is correct but in the wrong position; one number is incorrect and too large
    if not ((numbers[0] == '7' or numbers[1] == '7') and (numbers[0] != '8' and numbers[1] != '8')):
        return False
    # One letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet
    if not ((letters[0] == 'D' or letters[1] == 'D') and (letters[0] != 'I' and letters[1] != 'I')):
        return False

    # Check against the second guess: 50FH
    # Both numbers are incorrect and too small
    if not (int(numbers[0]) > 5 and int(numbers[1]) > 5):
        return False
    # One letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet
    if not ((letters[0] == 'H' or letters[1] == 'H') and (letters[0] != 'F' and letters[1] != 'F')):
        return False

    return True

# Generate all possible combinations of two numbers and two letters
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

for num_comb in permutations(numbers, 2):
    for let_comb in permutations(letters, 2):
        if is_valid_combination(num_comb, let_comb):
            password = list(num_comb) + list(let_comb)
            print(f"<<< {password} >>>")
            break