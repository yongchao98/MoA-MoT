from itertools import permutations

# Define possible domains
numbers = set(range(10))
letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Constraints based on feedback
def is_valid_combination(num1, num2, letter1, letter2):
    # Guess 1: 04XW
    if num1 in {0, 4} or num2 in {0, 4} or letter1 in "XW" or letter2 in "XW":
        return False
    # Guess 2: 06IH
    if num1 in {0, 6} or num2 in {0, 6} or (letter1 == 'I' and letter2 != 'H') or (letter2 == 'I' and letter1 != 'H'):
        return False
    # Guess 3: 54FA
    if num1 in {5, 4} or num2 in {5, 4} or letter1 in "FA" or letter2 in "FA":
        return False
    # Guess 4: 07ES
    if num1 in {0, 7} or num2 in {0, 7} or (letter1 == 'E' and letter2 != 'S') or (letter2 == 'E' and letter1 != 'S'):
        return False
    # Guess 5: 08HI
    if num1 != 8 or num2 in {0, 8} or letter1 != 'H' or letter2 in "HI":
        return False
    # Guess 6: 57OR
    if num1 in {5, 7} or num2 in {5, 7} or letter1 in "OR" or letter2 in "OR":
        return False
    # Guess 7: 29FB
    if num1 in {2, 9} or num2 != 2 or letter1 in "FB" or letter2 in "FB":
        return False
    return True

# Try all permutations of numbers and letters
for num1, num2 in permutations(numbers, 2):
    for letter1, letter2 in permutations(letters, 2):
        if is_valid_combination(num1, num2, letter1, letter2):
            password = [str(num1), str(num2), letter1, letter2]
            print(f"<<< {password} >>>")
            break