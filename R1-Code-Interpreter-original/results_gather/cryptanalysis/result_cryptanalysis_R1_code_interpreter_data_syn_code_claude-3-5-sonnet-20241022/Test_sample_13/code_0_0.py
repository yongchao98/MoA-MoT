def check_guess(actual, guess, feedback):
    # Unpack feedback
    num_correct_wrong_pos = feedback.get('num_correct_wrong_pos', 0)
    num_incorrect = feedback.get('num_incorrect', 0)
    letter_correct_pos = feedback.get('letter_correct_pos', 0)
    letter_too_early = feedback.get('letter_too_early', 0)
    letter_too_late = feedback.get('letter_too_late', 0)
    
    # Check numbers
    correct_wrong_pos = 0
    incorrect = 0
    for i in range(2):
        if guess[i] in actual[:2]:
            if guess[i] != actual[i]:
                correct_wrong_pos += 1
        else:
            incorrect += 1
            
    # Check letters
    correct_pos = 0
    too_early = 0
    too_late = 0
    for i in range(2):
        if guess[i+2] == actual[i+2]:
            correct_pos += 1
        elif guess[i+2] < actual[i+2]:
            too_early += 1
        elif guess[i+2] > actual[i+2]:
            too_late += 1
            
    return (correct_wrong_pos == num_correct_wrong_pos and
            incorrect == num_incorrect and
            correct_pos == letter_correct_pos and
            too_early == letter_too_early and
            too_late == letter_too_late)

def find_password():
    numbers = '0123456789'
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # Store all guesses and their feedback
    guesses = [
        ('95VB', {'num_correct_wrong_pos': 1, 'num_incorrect': 1, 'letter_correct_pos': 0}),
        ('34OH', {'num_correct_wrong_pos': 0, 'num_incorrect': 2, 'letter_correct_pos': 1, 'letter_too_early': 1}),
        ('63TK', {'num_correct_wrong_pos': 0, 'num_incorrect': 2, 'letter_correct_pos': 0}),
        ('61XN', {'num_correct_wrong_pos': 1, 'num_incorrect': 1, 'letter_correct_pos': 1, 'letter_too_late': 1})
    ]
    
    # Try all possible combinations
    for n1 in numbers:
        for n2 in numbers:
            if n1 != n2:  # Numbers can't repeat
                for l1 in letters:
                    for l2 in letters:
                        if l1 != l2:  # Letters can't repeat
                            candidate = n1 + n2 + l1 + l2
                            valid = True
                            
                            # Check against all guesses
                            for guess, feedback in guesses:
                                if not check_guess(candidate, guess, feedback):
                                    valid = False
                                    break
                            
                            if valid:
                                print([n1, n2, l1, l2])

find_password()