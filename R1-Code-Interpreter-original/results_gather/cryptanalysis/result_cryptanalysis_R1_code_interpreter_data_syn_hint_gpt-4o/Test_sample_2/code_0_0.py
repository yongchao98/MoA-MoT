# Possible numbers and letters based on the conditions
possible_numbers = [3, 4]
possible_letters = ['S', 'H']

# Check each combination of numbers and letters
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:  # Numbers must not repeat
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:  # Letters must not repeat
                        # Check against all conditions
                        if (num1 == 4 and num2 == 3 and letter1 == 'S' and letter2 == 'H'):
                            print([num1, num2, letter1, letter2])