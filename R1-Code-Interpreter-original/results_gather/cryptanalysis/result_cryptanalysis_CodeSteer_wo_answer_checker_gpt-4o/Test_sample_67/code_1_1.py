def is_valid_combination(numbers, letters, combination):
    num1, num2, letter1, letter2 = combination
    
    # Feedback rules
    feedbacks = [
        (79, 'FV', 'both numbers are incorrect and too large; both letters are incorrect'),
        (32, 'PZ', 'both numbers are incorrect and too small; both letters are incorrect and too late in the alphabet'),
        (9, 'EF', 'both numbers are incorrect; both letters are incorrect and too early in the alphabet'),
        (58, 'QD', 'one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect'),
        (79, 'OJ', 'both numbers are incorrect and too large; both letters are incorrect'),
        (64, 'LY', 'one number is correct but in the wrong position; one number is incorrect and too large; one letter is correct and in the correct position; one letter is incorrect and too late in the alphabet'),
        (48, 'HI', 'one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect and too early in the alphabet'),
        (45, 'TG', 'both numbers are correct and in the correct positions; both letters are incorrect'),
        (31, 'IB', 'both numbers are incorrect and too small; both letters are incorrect and too early in the alphabet'),
        (94, 'VW', 'one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect and too late in the alphabet'),
        (70, 'XN', 'both numbers are incorrect; both letters are incorrect and too late in the alphabet'),
        (70, 'BI', 'both numbers are incorrect; both letters are incorrect and too early in the alphabet'),
        (89, 'UG', 'both numbers are incorrect and too large; both letters are incorrect'),
        (70, 'KG', 'both numbers are incorrect; both letters are incorrect and too early in the alphabet')
    ]
    
    for guess, guess_letters, feedback in feedbacks:
        guess_num1, guess_num2 = divmod(guess, 10)
        
        if 'both numbers are incorrect and too large' in feedback:
            if num1 >= guess_num1 or num2 >= guess_num2:
                return False
        if 'both numbers are incorrect and too small' in feedback:
            if num1 <= guess_num1 or num2 <= guess_num2:
                return False
        if 'one number is correct but in the wrong position' in feedback:
            if not ((num1 == guess_num2 and num2 != guess_num1) or (num2 == guess_num1 and num1 != guess_num2)):
                return False
        if 'one number is correct and in the correct position' in feedback:
            if not ((num1 == guess_num1 and num2 != guess_num2) or (num2 == guess_num2 and num1 != guess_num1)):
                return False
        if 'both numbers are correct and in the correct positions' in feedback:
            if not (num1 == guess_num1 and num2 == guess_num2):
                return False
        
        if 'both letters are incorrect' in feedback:
            if letter1 in guess_letters or letter2 in guess_letters:
                return False
        if 'one letter is correct and in the correct position' in feedback:
            if not ((letter1 == guess_letters[0] and letter2 != guess_letters[1]) or (letter2 == guess_letters[1] and letter1 != guess_letters[0])):
                return False
        
        if 'both letters are incorrect and too early in the alphabet' in feedback:
            if letter1 < guess_letters[0] or letter2 < guess_letters[1]:
                return False
        if 'both letters are incorrect and too late in the alphabet' in feedback:
            if letter1 > guess_letters[0] or letter2 > guess_letters[1]:
                return False
        if 'one letter is incorrect and too late in the alphabet' in feedback:
            if letter1 > guess_letters[0] or letter2 > guess_letters[1]:
                return False
        if 'one letter is incorrect and too early in the alphabet' in feedback:
            if letter1 < guess_letters[0] or letter2 < guess_letters[1]:
                return False
    
    return True

def find_password():
    possible_numbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]
    
    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 == num2:
                continue
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 == letter2:
                        continue
                    combination = (num1, num2, letter1, letter2)
                    if is_valid_combination(possible_numbers, possible_letters, combination):
                        return combination

password = find_password()
print(f"<<< {list(map(str, password))} >>>")