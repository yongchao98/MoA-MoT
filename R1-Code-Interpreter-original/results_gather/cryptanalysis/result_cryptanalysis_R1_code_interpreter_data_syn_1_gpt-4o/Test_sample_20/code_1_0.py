# Possible numbers and letters based on deductions
possible_numbers = [0, 1, 3, 4]
possible_letters = ['A', 'B', 'D', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'X']

# Check each combination of numbers and letters
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:  # Numbers must not repeat
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:  # Letters must not repeat
                        # Check against all conditions
                        if (num1 == 3 or num2 == 3) and (num1 < 5 and num2 < 5) and \
                           (letter1 == 'X' or letter2 == 'X') and \
                           (letter1 < 'P' and letter2 < 'P'):
                            print([num1, num2, letter1, letter2])