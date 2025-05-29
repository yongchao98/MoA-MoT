from itertools import permutations

def check_guess(guess):
    # Check if a guess matches all conditions
    return all([
        check_14VM(guess),
        check_91UE(guess),
        check_85FO(guess),
        check_05QJ(guess),
        check_85ZQ(guess),
        check_95WT(guess),
        check_57KR(guess)
    ])

def check_14VM(guess):
    # One number correct and in position, one number too small, both letters wrong
    nums_correct = sum(1 for i in range(2) if guess[i] == '14'[i])
    return (nums_correct == 1 and 
            any(int(guess[i]) < int('14'[i]) for i in range(2)) and
            guess[2] not in 'VM' and guess[3] not in 'VM')

def check_91UE(guess):
    return (guess[0] != '9' and guess[1] != '1' and 
            guess[2] not in 'UE' and guess[3] not in 'UE')

def check_85FO(guess):
    return (guess[0] != '8' and guess[1] != '5' and 
            guess[2] not in 'FO' and guess[3] not in 'FO')

def check_05QJ(guess):
    return (guess[0] != '0' and guess[1] != '5' and
            ((guess[2] in 'QJ' and guess[2] != 'Q') or
             (guess[3] in 'QJ' and guess[3] != 'J')) and
            len(set(guess[2:]) & set('QJ')) == 1)

def check_85ZQ(guess):
    return (guess[0] != '8' and guess[1] != '5' and 
            guess[2] not in 'ZQ' and guess[3] not in 'ZQ')

def check_95WT(guess):
    return (guess[0] != '9' and guess[1] != '5' and 
            guess[2] not in 'WT' and guess[3] not in 'WT' and
            all(c < 'W' for c in guess[2:]))

def check_57KR(guess):
    return (guess[0] != '5' and guess[1] != '7' and 
            guess[2] not in 'KR' and guess[3] not in 'KR')

# Generate all possible combinations
# We can exclude numbers we know are wrong: 0,1,5,7,8,9
valid_numbers = [str(i) for i in range(10) if i not in {0,1,5,7,8,9}]
# We can exclude letters we know are wrong: V,M,U,E,F,O,Z,Q,W,T,K,R
invalid_letters = set('VMUEFOZQWTKR')
valid_letters = [chr(i) for i in range(65, 91) if chr(i) not in invalid_letters]

valid_solutions = []
for n1, n2 in permutations(valid_numbers, 2):
    for l1, l2 in permutations(valid_letters, 2):
        guess = n1 + n2 + l1 + l2
        if check_guess(guess):
            valid_solutions.append([n1, n2, l1, l2])

print(valid_solutions)