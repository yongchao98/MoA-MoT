from itertools import permutations

def is_valid_combination(combination):
    num1, num2, letter1, letter2 = combination
    
    # Debug: Print the current combination being checked
    print(f"Checking combination: {combination}")
    
    # Condition 1: 87VG
    # One number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect.
    if not ((num1 == '7' or num2 == '7') and (num1 != '8' and num2 != '8')):
        return False
    if letter1 in 'VG' or letter2 in 'VG':
        return False
    
    # Condition 2: 53KJ
    # Both numbers are incorrect; both letters are incorrect.
    if num1 in '53' or num2 in '53':
        return False
    if letter1 in 'KJ' or letter2 in 'KJ':
        return False
    
    # Condition 3: 68JZ
    # Both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and incorrect.
    if num1 in '68' or num2 in '68':
        return False
    if letter1 == 'Z' and letter2 != 'Z':
        return False
    if letter1 != 'Z' and letter2 == 'Z':
        return False
    if letter1 == 'J' or letter2 == 'J':
        return False
    
    # Condition 4: 46JT
    # Both numbers are incorrect; both letters are incorrect.
    if num1 in '46' or num2 in '46':
        return False
    if letter1 in 'JT' or letter2 in 'JT':
        return False
    
    # Condition 5: 39HZ
    # Both numbers are incorrect; both letters are correct and in the correct positions.
    if num1 in '39' or num2 in '39':
        return False
    if letter1 != 'H' or letter2 != 'Z':
        return False
    
    # Condition 6: 83UC
    # Both numbers are incorrect; both letters are incorrect.
    if num1 in '83' or num2 in '83':
        return False
    if letter1 in 'UC' or letter2 in 'UC':
        return False
    
    # Condition 7: 42CB
    # One number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect and too early in the alphabet.
    if not ((num1 == '2' and num2 != '4') or (num2 == '2' and num1 != '4')):
        return False
    if letter1 in 'CB' or letter2 in 'CB':
        return False
    
    return True

def find_password():
    numbers = '0123456789'
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    for num_comb in permutations(numbers, 2):
        for letter_comb in permutations(letters, 2):
            combination = num_comb + letter_comb
            if is_valid_combination(combination):
                return combination

password = find_password()
if password:
    print(f"<<< {list(password)} >>>")
else:
    print("No valid combination found.")