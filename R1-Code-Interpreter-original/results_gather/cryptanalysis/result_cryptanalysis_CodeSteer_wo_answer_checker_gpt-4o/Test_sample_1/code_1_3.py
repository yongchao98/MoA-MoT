def is_valid_combination(combination):
    num1, num2, letter1, letter2 = combination

    # Debug: Print the current combination being checked
    print(f"Checking combination: {combination}")

    # Condition 1: 92XL
    # One number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect.
    if not ((num1 == '9' or num2 == '9') and (num1 != '2' and num2 != '2') and letter1 not in 'XL' and letter2 not in 'XL'):
        print("Fails condition 1")
        return False

    # Condition 2: 75TE
    # Both numbers are incorrect; both letters are incorrect.
    if num1 in '75' or num2 in '75' or letter1 in 'TE' or letter2 in 'TE':
        print("Fails condition 2")
        return False

    # Condition 3: 38YF
    # One number is correct and in the correct position; one number is incorrect; both letters are incorrect.
    if not ((num1 == '3' and num2 != '8') or (num2 == '3' and num1 != '8')) or letter1 in 'YF' or letter2 in 'YF':
        print("Fails condition 3")
        return False

    # Condition 4: 14MX
    # Both numbers are incorrect; both letters are incorrect.
    if num1 in '14' or num2 in '14' or letter1 in 'MX' or letter2 in 'MX':
        print("Fails condition 4")
        return False

    # Condition 5: 56CY
    # Both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too late in the alphabet.
    if num1 in '56' or num2 in '56' or not (letter1 == 'C' or letter2 == 'C') or 'Y' in (letter1, letter2):
        print("Fails condition 5")
        return False

    # Condition 6: 09WV
    # One number is correct and in the correct position; one number is incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect.
    if not ((num1 == '0' and num2 != '9') or (num2 == '0' and num1 != '9')) or not ('W' in (letter1, letter2)) or 'V' in (letter1, letter2):
        print("Fails condition 6")
        return False

    return True

def find_password():
    possible_numbers = '0123456789'
    possible_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 == num2:
                continue
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 == letter2:
                        continue
                    combination = (num1, num2, letter1, letter2)
                    if is_valid_combination(combination):
                        return combination

    # If no valid combination is found, return an error message
    return "No valid combination found."

password = find_password()
print(f"<<< {list(password)} >>>")