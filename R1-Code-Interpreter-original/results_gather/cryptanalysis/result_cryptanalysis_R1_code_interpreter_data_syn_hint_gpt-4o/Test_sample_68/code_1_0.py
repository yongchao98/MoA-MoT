# Possible numbers and letters
correct_numbers = [1, 7]
possible_letters = ['C', 'D', 'E', 'F', 'I', 'J', 'O', 'Q', 'S', 'T', 'U', 'V', 'W', 'Z']

# Check combinations
for letter1 in possible_letters:
    for letter2 in possible_letters:
        if letter1 != letter2:
            # Check conditions
            if letter1 in ['C', 'D', 'E', 'F', 'I', 'J', 'O', 'Q', 'S', 'T', 'U', 'V', 'W', 'Z'] and letter2 in ['C', 'D', 'E', 'F', 'I', 'J', 'O', 'Q', 'S', 'T', 'U', 'V', 'W', 'Z']:
                print([1, 7, letter1, letter2])