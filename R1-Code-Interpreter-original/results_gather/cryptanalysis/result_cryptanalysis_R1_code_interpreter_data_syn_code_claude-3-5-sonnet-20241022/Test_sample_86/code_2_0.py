from itertools import permutations

def check_14VM(guess):
    # One number correct and in position
    # One number incorrect and too small
    n1, n2 = int(guess[0]), int(guess[1])
    correct_pos = (guess[0] == '1' and n2 != 4) or (guess[1] == '4' and n1 != 1)
    if guess[0] == '1':
        too_small = n2 < 4
    elif guess[1] == '4':
        too_small = n1 < 1
    else:
        return False
    return (correct_pos and too_small and 
            guess[2] not in 'VM' and guess[3] not in 'VM')

def check_91UE(guess):
    return (guess[0] != '9' and guess[1] != '1' and 
            guess[2] not in 'UE' and guess[3] not in 'UE')

def check_85FO(guess):
    return (guess[0] != '8' and guess[1] != '5' and 
            guess[2] not in 'FO' and guess[3] not in 'FO')

def check_05QJ(guess):
    # One letter correct but wrong position
    return (guess[0] != '0' and guess[1] != '5' and
            ((guess[2] in 'QJ' and guess[3] not in 'QJ') or
             (guess[3] in 'QJ' and guess[2] not in 'QJ')) and
            guess[2] != 'Q' and guess[3] != 'J')

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