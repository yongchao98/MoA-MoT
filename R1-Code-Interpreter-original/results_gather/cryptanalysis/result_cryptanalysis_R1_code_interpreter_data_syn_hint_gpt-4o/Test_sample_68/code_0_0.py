# Possible numbers and letters
possible_numbers = [1, 7]
possible_letters = ['C', 'D', 'E', 'F', 'I', 'J', 'O', 'Q', 'S', 'T', 'U', 'V', 'W', 'Z']

# Check combinations
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:
                        # Check conditions
                        if (num1 == 1 and num2 == 7) or (num1 == 7 and num2 == 1):
                            if letter1 in ['C', 'D', 'E', 'F', 'I', 'J', 'O', 'Q', 'S', 'T', 'U', 'V', 'W', 'Z'] and letter2 in ['C', 'D', 'E', 'F', 'I', 'J', 'O', 'Q', 'S', 'T', 'U', 'V', 'W', 'Z']:
                                print([num1, num2, letter1, letter2])