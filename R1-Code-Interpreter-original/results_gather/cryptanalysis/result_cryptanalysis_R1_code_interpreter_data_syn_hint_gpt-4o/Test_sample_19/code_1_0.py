# Possible numbers and letters based on the analysis
possible_numbers = [8, 9]
possible_letters = ['W', 'Z']

# Check the conditions
def check_combination(numbers, letters):
    # Condition 1: 18WV
    if numbers[0] == 8 and letters[0] == 'W':
        if numbers[1] > 1 and letters[1] > 'V':
            # Condition 2: 76GH
            if numbers[0] > 7 and numbers[1] > 6:
                if letters[0] > 'G' and letters[1] > 'H':
                    # Condition 3: 75LE
                    if numbers[0] > 7 and numbers[1] > 5:
                        if letters[0] > 'L' and letters[1] > 'E':
                            # Condition 4: 12BZ
                            if numbers[0] > 1 and numbers[1] > 2:
                                if letters[0] == 'Z' and letters[1] > 'B':
                                    return True
    return False

# Check all combinations
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:
                        if check_combination([num1, num2], [letter1, letter2]):
                            print([num1, num2, letter1, letter2])