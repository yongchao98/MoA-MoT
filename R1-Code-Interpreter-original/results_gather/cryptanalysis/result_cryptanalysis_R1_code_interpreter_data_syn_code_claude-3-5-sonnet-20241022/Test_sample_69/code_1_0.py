def check_guess(password, guess, feedback):
    # Unpack feedback
    num_correct_wrong_pos = feedback['num_correct_wrong_pos']
    num_too_large = feedback.get('num_too_large', 0)
    num_too_small = feedback.get('num_too_small', 0)
    num_correct_pos = feedback.get('num_correct_pos', 0)
    letter_correct_wrong_pos = feedback['letter_correct_wrong_pos']
    letter_too_late = feedback.get('letter_too_late', 0)
    letter_too_early = feedback.get('letter_too_early', 0)
    
    # Check numbers
    correct_pos = 0
    correct_wrong_pos = 0
    too_large = 0
    too_small = 0
    
    for i in range(2):
        if password[i] == guess[i]:
            correct_pos += 1
        elif guess[i] in password[:2]:
            correct_wrong_pos += 1
        elif int(guess[i]) > int(password[i]):
            too_large += 1
        else:
            too_small += 1
    
    # Check letters
    letter_correct_pos = 0
    letter_correct_wrong_pos = 0
    letter_late = 0
    letter_early = 0
    
    for i in range(2,4):
        if password[i] == guess[i]:
            letter_correct_pos += 1
        elif guess[i] in password[2:]:
            letter_correct_wrong_pos += 1
        elif guess[i] > password[i]:
            letter_late += 1
        else:
            letter_early += 1
            
    return (correct_pos == num_correct_pos and 
            correct_wrong_pos == num_correct_wrong_pos and
            too_large == num_too_large and
            too_small == num_too_small and
            letter_correct_wrong_pos == feedback['letter_correct_wrong_pos'] and
            letter_too_late == feedback.get('letter_too_late', 0) and
            letter_too_early == feedback.get('letter_too_early', 0))

# Test all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

feedbacks = [
    {'guess': '48TY', 
     'num_correct_wrong_pos': 1, 
     'num_too_large': 1,
     'letter_correct_wrong_pos': 1,
     'letter_too_late': 1},
    {'guess': '61LS',
     'num_correct_wrong_pos': 0,
     'num_too_large': 0,
     'letter_correct_wrong_pos': 0},
    {'guess': '21KP',
     'num_correct_pos': 1,
     'num_too_small': 1,
     'letter_correct_wrong_pos': 1,
     'letter_too_early': 1}
]

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # no repeating numbers
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # no repeating letters
                        password = [n1, n2, l1, l2]
                        valid = True
                        for feedback in feedbacks:
                            if not check_guess(password, feedback['guess'], feedback):
                                valid = False
                                break
                        if valid:
                            print(password)