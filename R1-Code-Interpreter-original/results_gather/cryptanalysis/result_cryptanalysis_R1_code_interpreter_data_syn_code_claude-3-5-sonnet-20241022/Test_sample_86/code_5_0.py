from itertools import permutations

def check_14VM(guess):
    return guess[0] == '1' and int(guess[1]) > 4

def check_91UE(guess):
    return (guess[0] != '9' and guess[1] != '1' and 
            guess[2] not in 'UE' and guess[3] not in 'UE')

def check_85FO(guess):
    return (guess[0] != '8' and guess[1] != '5' and 
            guess[2] not in 'FO' and guess[3] not in 'FO')

def check_05QJ(guess):
    return (guess[0] != '0' and guess[1] != '5' and
            guess[2] == 'J' and guess[3] not in 'QJ')  # J must be in first letter position

def check_85ZQ(guess):
    return (guess[0] != '8' and guess[1] != '5' and 
            guess[2] not in 'ZQ' and guess[3] not in 'ZQ')

def check_95WT(guess):
    return (guess[0] != '9' and guess[1] != '5' and 
            all(c < 'W' for c in guess[2:]))

def check_57KR(guess):
    return (guess[0] != '5' and guess[1] != '7' and 
            guess[2] not in 'KR' and guess[3] not in 'KR')

def check_letter_order(guess):
    # Second letter must be after J in alphabet
    return guess[3] > 'J'

# We know first number is 1 and second is 6
numbers = ['1', '6']
# We know first letter is J
# Second letter must be after J but before W, and not in invalid set
invalid_letters = set('VMUEFOZKRWTQ')
valid_second_letters = [c for c in [chr(i) for i in range(ord('J')+1, ord('W'))]
                       if c not in invalid_letters]

solutions = []
for l2 in valid_second_letters:
    guess = numbers[0] + numbers[1] + 'J' + l2
    if all([
        check_14VM(guess),
        check_91UE(guess),
        check_85FO(guess),
        check_05QJ(guess),
        check_85ZQ(guess),
        check_95WT(guess),
        check_57KR(guess),
        check_letter_order(guess)
    ]):
        solutions.append([numbers[0], numbers[1], 'J', l2])

print(solutions)