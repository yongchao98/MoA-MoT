def is_valid_combination(numbers, letters, combination):
    num1, num2, letter1, letter2 = combination
    
    # Condition 1: 25SN
    if num1 in [2, 5] or num2 in [2, 5] or letter1 in ['S', 'N'] or letter2 in ['S', 'N']:
        return False
    
    # Condition 2: 28CK
    if num1 in [2, 8] or num2 in [2, 8] or letter1 in ['C', 'K'] or letter2 in ['C', 'K']:
        return False
    
    # Condition 3: 93CH
    if not (num1 == 9 and num2 != 3) or letter1 in ['C', 'H'] or letter2 in ['C', 'H']:
        return False
    
    # Condition 4: 24RD
    if num1 in [2, 4] or num2 in [2, 4] or (letter1 == 'R' and letter2 != 'D') or letter2 == 'R':
        return False
    
    # Condition 5: 07DM
    if not (num2 == 0 and num1 != 7) or letter1 != 'D' or letter2 != 'M':
        return False
    
    return True

# Possible numbers and letters based on feedback
possible_numbers = [0, 1, 3, 6, 9]
possible_letters = ['A', 'B', 'D', 'E', 'M']

# Iterate through all possible combinations
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
                    print(f"<<< {list(combination)} >>>")