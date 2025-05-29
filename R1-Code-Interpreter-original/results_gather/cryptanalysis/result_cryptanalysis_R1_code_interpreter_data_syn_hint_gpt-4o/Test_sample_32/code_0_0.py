# Possible numbers and letters based on the analysis
possible_numbers = [0, 8]
possible_letters = ['A', 'E']

# Check all combinations of the possible numbers and letters
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:  # Numbers must not repeat
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:  # Letters must not repeat
                        # Check against all conditions
                        if (num1 == 8 or num2 == 8) and (letter1 == 'A' or letter2 == 'A'):
                            print([num1, num2, letter1, letter2])