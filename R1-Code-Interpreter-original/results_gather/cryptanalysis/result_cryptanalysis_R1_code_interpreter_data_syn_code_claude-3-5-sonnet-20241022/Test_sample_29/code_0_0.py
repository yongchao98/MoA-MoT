def check_guess(guess, actual):
    feedback = {'correct_pos_num': 0, 'wrong_pos_num': 0, 'small_num': 0,
                'correct_pos_letter': 0, 'wrong_pos_letter': 0, 'early_letter': 0}
    
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            feedback['correct_pos_num'] += 1
        elif guess[i] in actual[:2]:
            feedback['wrong_pos_num'] += 1
        elif int(guess[i]) < int(actual[i]):
            feedback['small_num'] += 1
            
    # Check letters
    for i in range(2,4):
        if guess[i] == actual[i]:
            feedback['correct_pos_letter'] += 1
        elif guess[i] in actual[2:]:
            feedback['wrong_pos_letter'] += 1
        elif guess[i] < actual[i]:
            feedback['early_letter'] += 1
            
    return feedback

def matches_conditions(candidate):
    # Guess 1: 74JY
    f1 = check_guess('74JY', candidate)
    if not (f1['correct_pos_num'] == 1 and f1['small_num'] == 1 and 
            f1['wrong_pos_letter'] == 1 and f1['early_letter'] == 1):
        return False
    
    # Guess 2: 93ZN
    f2 = check_guess('93ZN', candidate)
    if not (f2['correct_pos_num'] == 0 and f2['wrong_pos_num'] == 0 and 
            f2['correct_pos_letter'] == 0 and f2['wrong_pos_letter'] == 0):
        return False
    
    # Guess 3: 26MU
    f3 = check_guess('26MU', candidate)
    if not (f3['small_num'] == 2 and 
            f3['correct_pos_letter'] == 0 and f3['wrong_pos_letter'] == 0):
        return False
    
    # Guess 4: 57FS
    f4 = check_guess('57FS', candidate)
    if not (f4['wrong_pos_num'] == 1 and f4['small_num'] == 1 and 
            f4['correct_pos_letter'] == 1 and f4['early_letter'] == 1):
        return False
    
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
found = []

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        candidate = n1 + n2 + l1 + l2
                        if matches_conditions(candidate):
                            found.append(candidate)

print(found)