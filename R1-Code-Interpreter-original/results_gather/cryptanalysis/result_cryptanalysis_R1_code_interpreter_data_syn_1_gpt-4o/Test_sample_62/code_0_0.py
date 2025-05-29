# Possible numbers and letters
possible_numbers = [0, 6, 8]
possible_letters = ['E', 'F', 'H', 'I', 'N', 'O', 'P', 'Q', 'R', 'W']

# We know 6 is in the correct position from guess 6
# We need to find the other number and the two letters

# Iterate over possible numbers and letters to find the correct combination
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:  # Numbers must not repeat
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:  # Letters must not repeat
                        # Check conditions
                        if (num1 == 6 or num2 == 6) and (num1 != 6 or num2 != 6):
                            if (letter1 == 'I' or letter2 == 'I') and (letter1 != 'W' and letter2 != 'W'):
                                print([num1, num2, letter1, letter2])