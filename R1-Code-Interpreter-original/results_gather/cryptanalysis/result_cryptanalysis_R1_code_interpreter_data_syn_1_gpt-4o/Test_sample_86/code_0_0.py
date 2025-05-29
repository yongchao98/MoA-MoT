# Possible numbers and letters based on the deductions
possible_numbers = [4, 6]
possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

# Iterate over possible combinations
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:
                        # Check conditions
                        if (num1 == 4 and letter2 == 'J'):
                            print([num1, num2, letter1, letter2])