from itertools import permutations

def is_valid_combination(combo, guesses):
    for guess, feedback in guesses:
        numbers_match = (combo[0] == guess[0] and combo[1] == guess[1])
        letters_match = (combo[2] == guess[2] and combo[3] == guess[3])
        
        if feedback == "both numbers are correct and in the correct positions":
            if not numbers_match:
                return False
        elif feedback == "both numbers are incorrect":
            if numbers_match:
                return False
        elif feedback == "one number is correct and in the correct position":
            if not ((combo[0] == guess[0] or combo[1] == guess[1]) and not numbers_match):
                return False
        
        if feedback == "one letter is correct and in the correct position":
            if not (combo[2] == guess[2] or combo[3] == guess[3]):
                return False
        elif feedback == "both letters are incorrect":
            if letters_match:
                return False
        elif feedback == "one letter is correct but in the wrong position":
            if not ((combo[2] == guess[3] or combo[3] == guess[2]) and not letters_match):
                return False

    return True

def deduce_password():
    numbers = ['2', '8']
    letters = ['M', 'V', 'I', 'P', 'T']
    guesses = [
        ("28MV", "both numbers are correct and in the correct positions; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet"),
        ("41FX", "both numbers are incorrect; both letters are incorrect"),
        ("58VP", "one number is correct and in the correct position; one number is incorrect and incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet"),
        ("93IT", "both numbers are incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet")
    ]

    for num_combo in permutations(numbers, 2):
        for letter_combo in permutations(letters, 2):
            combo = num_combo + letter_combo
            if is_valid_combination(combo, guesses):
                return combo

password = deduce_password()
print(f"<<< {list(password)} >>>")