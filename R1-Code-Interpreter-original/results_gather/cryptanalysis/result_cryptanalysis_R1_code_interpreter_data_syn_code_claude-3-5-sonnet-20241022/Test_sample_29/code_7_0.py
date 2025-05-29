def check_guess(candidate, guess, expected):
    """
    expected is a tuple of:
    (correct_pos_nums, wrong_pos_nums, small_nums, 
     correct_pos_letters, wrong_pos_letters, early_letters)
    """
    result = [0, 0, 0, 0, 0, 0]
    
    # Check numbers (positions 0 and 1)
    for i in range(2):
        if guess[i] == candidate[i]:
            result[0] += 1  # correct position number
        elif guess[i] in candidate[:2]:
            result[1] += 1  # wrong position number
        elif int(guess[i]) < int(candidate[i]):
            result[2] += 1  # too small number
            
    # Check letters (positions 2 and 3)
    for i in range(2, 4):
        if guess[i] == candidate[i]:
            result[3] += 1  # correct position letter
        elif guess[i] in candidate[2:]:
            result[4] += 1  # wrong position letter
        elif guess[i] < candidate[i]:
            result[5] += 1  # too early letter
            
    return tuple(result) == expected

def find_password():
    numbers = '0123456789'
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    candidates = []
    
    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    
                    pwd = n1 + n2 + l1 + l2
                    
                    # Check all conditions
                    
                    # Guess 1: 74JY
                    c1 = check_guess(pwd, "74JY", (1, 0, 1, 0, 1, 1))
                    
                    # Guess 2: 93ZN
                    c2 = check_guess(pwd, "93ZN", (0, 0, 0, 0, 0, 0))
                    
                    # Guess 3: 26MU
                    c3 = check_guess(pwd, "26MU", (0, 0, 2, 0, 0, 0))
                    
                    # Guess 4: 57FS
                    c4 = check_guess(pwd, "57FS", (0, 1, 1, 1, 0, 1))
                    
                    # Additional checks
                    c5 = ('9' not in pwd and '3' not in pwd and 
                          'Z' not in pwd and 'N' not in pwd)
                    c6 = ('2' not in pwd and '6' not in pwd and 
                          'M' not in pwd and 'U' not in pwd)
                    
                    if all([c1, c2, c3, c4, c5, c6]):
                        candidates.append(pwd)
    
    return candidates

result = find_password()
print(result)