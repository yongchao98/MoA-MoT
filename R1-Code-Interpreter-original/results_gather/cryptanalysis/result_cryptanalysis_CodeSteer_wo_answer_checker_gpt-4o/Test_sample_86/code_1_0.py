def is_valid_combination(numbers, letters, combination):
    num1, num2, letter1, letter2 = combination

    # Check each condition
    conditions = [
        (lambda: (num1 == '1' and num2 != '4') or (num1 != '1' and num2 == '4'), '14VM'),
        (lambda: num1 != '9' and num2 != '1', '91UE'),
        (lambda: num1 != '8' and num2 != '5', '85FO'),
        (lambda: num1 != '0' and num2 != '5' and (letter1 == 'J' or letter2 == 'J'), '05QJ'),
        (lambda: num1 != '8' and num2 != '5', '85ZQ'),
        (lambda: num1 != '9' and num2 != '5' and letter1 < 'W' and letter2 < 'T', '95WT'),
        (lambda: num1 != '5' and num2 != '7', '57KR')
    ]

    for condition, guess in conditions:
        if not condition():
            return False

    return True

def find_password():
    possible_numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

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
print(f'<<< {list(password)} >>>')