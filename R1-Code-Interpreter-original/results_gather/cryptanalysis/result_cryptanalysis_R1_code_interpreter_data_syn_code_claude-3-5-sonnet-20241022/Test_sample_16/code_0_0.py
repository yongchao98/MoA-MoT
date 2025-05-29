def check_guess(actual, guess, feedback):
    # Extract numbers and letters
    actual_nums = [actual[0], actual[1]]
    actual_lets = [actual[2], actual[3]]
    guess_nums = [guess[0], guess[1]]
    guess_lets = [guess[2], guess[3]]
    
    # Check numbers
    if feedback['nums'] == 'both_incorrect_small':
        if not (int(guess_nums[0]) < int(actual_nums[0]) and int(guess_nums[1]) < int(actual_nums[1])):
            return False
    elif feedback['nums'] == 'both_incorrect':
        if guess_nums[0] in actual_nums or guess_nums[1] in actual_nums:
            return False
    elif feedback['nums'] == 'one_correct_wrong_pos_one_small':
        found = False
        for i in range(2):
            for j in range(2):
                if i != j and guess_nums[i] == actual_nums[j] and int(guess_nums[1-i]) < int(actual_nums[1-i]):
                    found = True
        if not found:
            return False
    elif feedback['nums'] == 'one_correct_wrong_pos_one_incorrect':
        found = False
        for i in range(2):
            for j in range(2):
                if i != j and guess_nums[i] == actual_nums[j] and guess_nums[1-i] != actual_nums[0] and guess_nums[1-i] != actual_nums[1]:
                    found = True
        if not found:
            return False
    
    # Check letters
    if feedback['lets'] == 'both_incorrect':
        if guess_lets[0] in actual_lets or guess_lets[1] in actual_lets:
            return False
    elif feedback['lets'] == 'one_correct_pos_one_late':
        correct = False
        late = False
        for i in range(2):
            if guess_lets[i] == actual_lets[i]:
                correct = True
            if ord(guess_lets[i]) > ord(actual_lets[i]):
                late = True
        if not (correct and late):
            return False
    elif feedback['lets'] == 'one_correct_pos':
        correct = sum(1 for i in range(2) if guess_lets[i] == actual_lets[i])
        if correct != 1:
            return False
    
    return True

def find_password():
    numbers = [str(i).zfill(2) for i in range(10)]
    letters = [chr(i) for i in range(65, 91)]  # A-Z
    
    # Test all possible combinations
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    
                    candidate = [str(n1), str(n2), l1, l2]
                    
                    # Check against all guesses
                    if not check_guess(candidate, ['8', '0', 'K', 'R'], 
                                    {'nums': 'one_correct_wrong_pos_one_small', 'lets': 'both_incorrect'}):
                        continue
                    
                    if not check_guess(candidate, ['8', '7', 'N', 'P'],
                                    {'nums': 'one_correct_wrong_pos_one_incorrect', 'lets': 'one_correct_pos'}):
                        continue
                    
                    if not check_guess(candidate, ['4', '7', 'A', 'Z'],
                                    {'nums': 'both_incorrect', 'lets': 'one_correct_pos_one_late'}):
                        continue
                    
                    if not check_guess(candidate, ['1', '5', 'R', 'F'],
                                    {'nums': 'both_incorrect_small', 'lets': 'both_incorrect'}):
                        continue
                    
                    print(candidate)

find_password()