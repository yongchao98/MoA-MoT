def is_valid_combination(numbers, letters, feedback):
    for guess, response in feedback:
        guess_numbers, guess_letters = guess[:2], guess[2:]
        
        # Check numbers
        if response['numbers'] == 'both incorrect':
            if any(num in numbers for num in guess_numbers):
                return False
        elif response['numbers'] == 'both correct and in correct positions':
            if numbers != guess_numbers:
                return False
        
        # Check letters
        if response['letters'] == 'one correct and in correct position':
            if not any(letters[i] == guess_letters[i] for i in range(2)):
                return False
        elif response['letters'] == 'one correct but in wrong position':
            if not any(letters[i] != guess_letters[i] and letters[i] in guess_letters for i in range(2)):
                return False
        elif response['letters'] == 'both incorrect':
            if any(letter in letters for letter in guess_letters):
                return False
        elif response['letters'] == 'both incorrect and too early in the alphabet':
            if any(letter in letters or letter <= guess_letters[i] for i, letter in enumerate(letters)):
                return False

    return True

def find_password():
    feedback = [
        ("18VC", {'numbers': 'both incorrect', 'letters': 'one correct and in correct position'}),
        ("54OF", {'numbers': 'both incorrect', 'letters': 'both incorrect and too early in the alphabet'}),
        ("36AV", {'numbers': 'both correct and in correct positions', 'letters': 'one correct but in wrong position'}),
        ("42GE", {'numbers': 'both incorrect', 'letters': 'both incorrect and too early in the alphabet'}),
        ("08NQ", {'numbers': 'both incorrect', 'letters': 'both incorrect'})
    ]

    # Possible numbers and letters based on feedback
    possible_numbers = ['3', '6']
    possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1) if i > ord('V')]

    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 == num2:
                continue
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 == letter2:
                        continue
                    if is_valid_combination([num1, num2], [letter1, letter2], feedback):
                        return [num1, num2, letter1, letter2]

password = find_password()
print(f"<<< {password} >>>")