def is_valid_combination(numbers, letters, combination):
    # Unpack the combination
    num1, num2, letter1, letter2 = combination
    
    # Check each condition
    # Condition 1: 47YQ
    if num1 >= 4 or num2 >= 4 or letter1 >= 'Y' or letter2 >= 'Q':
        return False
    
    # Condition 2: 91LQ
    if num1 in [9, 1] or num2 in [9, 1] or (letter1 != 'L' and letter2 != 'L') or letter2 >= 'Q':
        return False
    
    # Condition 3: 23TP
    if (num1 != 2 and num2 != 2) or (num1 == 3 or num2 == 3) or letter1 >= 'T' or letter2 >= 'P':
        return False
    
    # Condition 4: 09KG
    if (num1 != 0 and num2 != 0) or num1 >= 9 or (letter1 != 'K' and letter2 != 'K') or letter2 == 'G':
        return False
    
    return True

def find_password():
    possible_numbers = [0, 2, 3]  # Based on feedback
    possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']  # Letters before L
    
    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 == num2:
                continue
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 == letter2:
                        continue
                    combination = (num1, num2, letter1, letter2)
                    if is_valid_combination(possible_numbers, possible_letters, combination):
                        return combination

password = find_password()
print(f"<<< {list(password)} >>>")