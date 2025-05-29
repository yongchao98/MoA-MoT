def check_guess(guess, actual):
    # Convert both to lists for easier comparison
    guess = list(guess)
    actual = list(actual)
    
    # Check conditions for each guess
    def check_41CK(test):
        correct_pos_num = 0
        correct_wrong_pos_num = 0
        correct_pos_letter = 0
        correct_wrong_pos_letter = 0
        small_num = 0
        early_letter = 0
        
        # Check numbers
        for i in [0, 1]:
            if guess[i] == test[i]:
                correct_pos_num += 1
            elif int(guess[i]) < int(test[i]):
                small_num += 1
                
        # Check letters
        for i in [2, 3]:
            if guess[i] == test[i]:
                correct_pos_letter += 1
            elif guess[i] < test[i]:
                early_letter += 1
            elif guess[i] in test:
                correct_wrong_pos_letter += 1
                
        return (correct_pos_num == 1 and small_num == 1 and 
                correct_wrong_pos_letter == 1 and 
                (correct_pos_letter + correct_wrong_pos_letter + early_letter == 2))

    def check_64DE(test):
        correct_wrong_pos_num = 0
        for i in [0, 1]:
            if guess[i] in test[:2] and guess[i] != test[i]:
                correct_wrong_pos_num += 1
        
        early_letters = sum(1 for i in [2, 3] if guess[i] < test[i])
        return correct_wrong_pos_num == 1 and early_letters == 2

    def check_87JY(test):
        return not any(n in test for n in ['8', '7']) and not any(l in test for l in ['J', 'Y'])

    def check_12OD(test):
        small_nums = sum(1 for i in [0, 1] if int(guess[i]) < int(test[i]))
        correct_wrong_pos_letter = sum(1 for i in [2, 3] if guess[i] in test[2:] and guess[i] != test[i])
        early_letter = sum(1 for i in [2, 3] if guess[i] < test[i] and guess[i] not in test[2:])
        return small_nums == 2 and correct_wrong_pos_letter == 1 and early_letter == 1

    if guess == '41CK':
        return check_41CK(actual)
    elif guess == '64DE':
        return check_64DE(actual)
    elif guess == '87JY':
        return check_87JY(actual)
    elif guess == '12OD':
        return check_12OD(actual)
    return False

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
solutions = []

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        test = n1 + n2 + l1 + l2
                        if (check_guess('41CK', test) and 
                            check_guess('64DE', test) and 
                            check_guess('87JY', test) and 
                            check_guess('12OD', test)):
                            solutions.append([n1, n2, l1, l2])

print(solutions)