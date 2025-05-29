def is_valid_guess(candidate, guess, feedback):
    # Unpack feedback
    num_correct_pos = sum(1 for i in range(2) if feedback[i] == 'correct_pos' and candidate[i] == guess[i])
    num_too_large = sum(1 for i in range(2) if feedback[i] == 'too_large' and int(guess[i]) > int(candidate[i]))
    letter_too_late = sum(1 for i in range(2,4) if feedback[i] == 'too_late' and guess[i] > candidate[i])
    letter_correct_pos = sum(1 for i in range(2,4) if feedback[i] == 'correct_pos' and candidate[i] == guess[i])
    letter_wrong = sum(1 for i in range(2,4) if feedback[i] == 'wrong' and candidate[i] != guess[i])
    num_wrong_pos = sum(1 for i in range(2) if feedback[i] == 'wrong_pos' and guess[i] in candidate[:2] and guess[i] != candidate[i])
    letter_wrong_pos = sum(1 for i in range(2,4) if feedback[i] == 'wrong_pos' and guess[i] in candidate[2:] and guess[i] != candidate[i])
    
    return all([
        num_correct_pos == sum(1 for i in range(2) if feedback[i] == 'correct_pos'),
        num_too_large == sum(1 for i in range(2) if feedback[i] == 'too_large'),
        letter_too_late == sum(1 for i in range(2,4) if feedback[i] == 'too_late'),
        letter_correct_pos == sum(1 for i in range(2,4) if feedback[i] == 'correct_pos'),
        letter_wrong == sum(1 for i in range(2,4) if feedback[i] == 'wrong'),
        num_wrong_pos == sum(1 for i in range(2) if feedback[i] == 'wrong_pos'),
        letter_wrong_pos == sum(1 for i in range(2,4) if feedback[i] == 'wrong_pos')
    ])

def find_password():
    guesses = [
        ('08VY', ['correct_pos', 'too_large', 'too_late', 'too_late']),
        ('07CX', ['correct_pos', 'too_large', 'wrong', 'wrong']),
        ('45NJ', ['too_large', 'too_large', 'correct_pos', 'wrong']),
        ('28YQ', ['wrong_pos', 'too_large', 'wrong_pos', 'too_late'])
    ]
    
    import itertools
    numbers = '0123456789'
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # Try all possible combinations
    for num_pair in itertools.permutations(numbers, 2):
        for letter_pair in itertools.permutations(letters, 2):
            candidate = list(num_pair + letter_pair)
            if all(is_valid_guess(candidate, list(guess[0]), guess[1]) for guess in guesses):
                return candidate
    
    return None

result = find_password()
print(result)