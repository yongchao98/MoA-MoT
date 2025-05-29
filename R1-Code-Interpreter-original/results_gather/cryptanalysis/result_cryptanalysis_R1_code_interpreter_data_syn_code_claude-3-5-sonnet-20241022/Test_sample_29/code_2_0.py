def check_guess(guess, actual):
    # Split into numbers and letters
    guess_nums = [guess[0], guess[1]]
    guess_letters = [guess[2], guess[3]]
    actual_nums = [actual[0], actual[1]]
    actual_letters = [actual[2], actual[3]]
    
    result = {
        'correct_pos_num': 0,
        'wrong_pos_num': 0,
        'small_num': 0,
        'correct_pos_letter': 0,
        'wrong_pos_letter': 0,
        'early_letter': 0
    }
    
    # Check numbers
    for i in range(2):
        if guess_nums[i] == actual_nums[i]:
            result['correct_pos_num'] += 1
        elif guess_nums[i] in actual_nums:
            result['wrong_pos_num'] += 1
        elif int(guess_nums[i]) < int(actual_nums[i]):
            result['small_num'] += 1
    
    # Check letters
    for i in range(2):
        if guess_letters[i] == actual_letters[i]:
            result['correct_pos_letter'] += 1
        elif guess_letters[i] in actual_letters:
            result['wrong_pos_letter'] += 1
        elif guess_letters[i] < actual_letters[i]:
            result['early_letter'] += 1
    
    return result

def matches_all_conditions(password):
    # Guess 1: 74JY
    f1 = check_guess('74JY', password)
    c1 = (f1['correct_pos_num'] == 1 and 
          f1['small_num'] == 1 and 
          f1['wrong_pos_letter'] == 1 and 
          f1['early_letter'] == 1)
    
    # Guess 2: 93ZN
    f2 = check_guess('93ZN', password)
    c2 = (f2['correct_pos_num'] == 0 and 
          f2['wrong_pos_num'] == 0 and 
          f2['correct_pos_letter'] == 0 and 
          f2['wrong_pos_letter'] == 0)
    
    # Guess 3: 26MU
    f3 = check_guess('26MU', password)
    c3 = (f3['small_num'] == 2 and 
          f3['correct_pos_letter'] == 0 and 
          f3['wrong_pos_letter'] == 0)
    
    # Guess 4: 57FS
    f4 = check_guess('57FS', password)
    c4 = (f4['wrong_pos_num'] == 1 and 
          f4['small_num'] == 1 and 
          f4['correct_pos_letter'] == 1 and 
          f4['early_letter'] == 1)
    
    return c1 and c2 and c3 and c4

# Generate and test all possible combinations
valid_passwords = []
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        password = n1 + n2 + l1 + l2
                        if matches_all_conditions(password):
                            valid_passwords.append(password)

print(valid_passwords)