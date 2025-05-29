from itertools import permutations

def is_valid_combination(numbers, letters, guess, feedback):
    guess_numbers, guess_letters = guess[:2], guess[2:]
    feedback_numbers, feedback_letters = feedback[:2], feedback[2:]
    
    # Check numbers
    for i, (g, f) in enumerate(zip(guess_numbers, feedback_numbers)):
        if f == 'correct':
            if numbers[i] != g:
                return False
        elif f == 'wrong position':
            if numbers[i] == g or g not in numbers:
                return False
        elif f == 'incorrect':
            if g in numbers:
                return False
    
    # Check letters
    for i, (g, f) in enumerate(zip(guess_letters, feedback_letters)):
        if f == 'correct':
            if letters[i] != g:
                return False
        elif f == 'wrong position':
            if letters[i] == g or g not in letters:
                return False
        elif f == 'incorrect':
            if g in letters:
                return False
    
    return True

def deduce_password():
    possible_numbers = [str(i) for i in range(10)]
    possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]
    
    # Feedback for each guess
    feedbacks = [
        (("8", "6", "W", "U"), ("incorrect", "incorrect", "wrong position", "incorrect")),
        (("9", "4", "M", "F"), ("wrong position", "incorrect", "incorrect", "incorrect")),
        (("8", "0", "U", "V"), ("incorrect", "incorrect", "incorrect", "incorrect")),
        (("4", "7", "T", "H"), ("incorrect", "incorrect", "correct", "incorrect")),
        (("0", "3", "S", "U"), ("wrong position", "incorrect", "incorrect", "incorrect"))
    ]
    
    # Generate all possible combinations of two numbers and two letters
    for numbers in permutations(possible_numbers, 2):
        for letters in permutations(possible_letters, 2):
            if all(is_valid_combination(numbers, letters, guess, feedback) for guess, feedback in feedbacks):
                print(f"Valid combination found: {numbers + letters}")
                return list(numbers) + list(letters)

password = deduce_password()
print(f"<<< {password} >>>")