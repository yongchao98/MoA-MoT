# Possible numbers and letters based on deductions
numbers = [3, 4]
letters = ['S', 'F', 'C', 'P']

# Check each combination of numbers and letters
for num1 in numbers:
    for num2 in numbers:
        if num1 != num2:  # Numbers must not repeat
            for letter1 in letters:
                for letter2 in letters:
                    if letter1 != letter2:  # Letters must not repeat
                        # Check conditions
                        if (num1 == 3 and num2 == 4) and (letter1 == 'S' or letter2 == 'S') and (letter1 == 'C' or letter2 == 'C'):
                            print([num1, num2, letter1, letter2])