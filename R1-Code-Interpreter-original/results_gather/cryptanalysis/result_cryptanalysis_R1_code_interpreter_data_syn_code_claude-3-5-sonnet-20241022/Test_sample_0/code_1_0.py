def check_guess(password, guess, feedback):
    # Unpack feedback
    num_correct_wrong_pos = sum(1 for i in range(2) if 'one number is correct but in the wrong position' in feedback and 
                               str(guess[i]) in str(password) and str(guess[i]) != str(password[i]))
    num_too_large = sum(1 for i in range(2) if 'too large' in feedback and 
                       int(guess[i]) > int(password[i]))
    letter_correct_pos = sum(1 for i in range(2,4) if 'correct and in the correct position' in feedback and 
                           guess[i] == password[i])
    letter_correct_wrong_pos = sum(1 for i in range(2,4) if 'correct but in the wrong position' in feedback and 
                                 guess[i] in password[2:] and guess[i] != password[i])
    letter_too_late = sum(1 for i in range(2,4) if 'too late' in feedback and 
                         guess[i] > password[i])
    
    return (
        ('one number is correct but in the wrong position' in feedback) == (num_correct_wrong_pos > 0) and
        ('too large' in feedback) == (num_too_large > 0) and
        ('correct and in the correct position' in feedback) == (letter_correct_pos > 0) and
        ('correct but in the wrong position' in feedback) == (letter_correct_wrong_pos > 0) and
        ('too late' in feedback) == (letter_too_late > 0)
    )

def find_password():
    # Generate all possible numbers and letters
    numbers = list(range(0, 10))
    letters = [chr(i) for i in range(65, 91)]  # A-Z
    
    # Test cases
    guesses = [
        ('38TA', 'one number is correct but in the wrong position; one number is incorrect and too large; one letter is correct and in the correct position; one letter is incorrect and too late in the alphabet'),
        ('10GY', 'both numbers are incorrect and too small; both letters are incorrect and too late in the alphabet'),
        ('19QH', 'both numbers are incorrect; both letters are incorrect and too late in the alphabet'),
        ('12EF', 'both numbers are incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect and incorrect'),
        ('65TL', 'one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect and too late in the alphabet')
    ]
    
    # Try all possible combinations
    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:  # Skip if numbers repeat
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:  # Skip if letters repeat
                        continue
                    password = [str(n1), str(n2), l1, l2]
                    
                    # Check if this password satisfies all guesses
                    if all(check_guess(password, list(guess[0]), guess[1]) for guess in guesses):
                        print(password)
                        return

find_password()