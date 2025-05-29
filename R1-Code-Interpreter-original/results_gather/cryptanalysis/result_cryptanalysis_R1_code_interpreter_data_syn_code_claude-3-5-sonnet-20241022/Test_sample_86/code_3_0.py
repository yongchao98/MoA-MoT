from itertools import permutations

def check_14VM(guess):
    # One number correct and in position
    # One number incorrect and too small
    n1, n2 = int(guess[0]), int(guess[1])
    if guess[0] == '1':  # If 1 is correct
        return n2 > 4 and guess[2] not in 'VM' and guess[3] not in 'VM'
    elif guess[1] == '4':  # If 4 is correct
        return n1 < 1 and guess[2] not in 'VM' and guess[3] not in 'VM'
    return False

def check_91UE(guess):
    return (guess[0] != '9' and guess[1] != '1' and 
            guess[2] not in 'UE' and guess[3] not in 'UE')

def check_85FO(guess):
    return (guess[0] != '8' and guess[1] != '5' and 
            guess[2] not in 'FO' and guess[3] not in 'FO')

def check_05QJ(guess):
    # One letter correct but wrong position
    letters = set(guess[2:])
    return (guess[0] != '0' and guess[1] != '5' and
            'J' in letters and  # J must be one of the letters
            guess[2] != 'Q' and guess[3] != 'J' and  # but in different position
            len(letters & set('QJ')) == 1)  # only one of Q or J should be present

def check_85ZQ(guess):
    return (guess[0] != '8' and guess[1] != '5' and 
            guess[2] not in 'ZQ' and guess[3] not in 'ZQ')

def check_95WT(guess):
    return (guess[0] != '9' and guess[1] != '5' and 
            all(c < 'W' for c in guess[2:]))

def check_57KR(guess):
    return (guess[0] != '5' and guess[1] != '7' and 
            guess[2] not in 'KR' and guess[3] not in 'KR')

def check_all(guess):
    return all([
        check_14VM(guess),
        check_91UE(guess),
        check_85FO(guess),
        check_05QJ(guess),
        check_85ZQ(guess),
        check_95WT(guess),
        check_57KR(guess)
    ])

# Try all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z

solutions = []
for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(letters, 2):
        guess = n1 + n2 + l1 + l2
        if check_all(guess):
            solutions.append([n1, n2, l1, l2])

print(solutions)