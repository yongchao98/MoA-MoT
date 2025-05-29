from itertools import permutations

def check_guess(guess, actual):
    # Check if a guess matches all the feedback conditions
    numbers = [guess[0], guess[1]]
    letters = [guess[2], guess[3]]
    actual_numbers = [actual[0], actual[1]]
    actual_letters = [actual[2], actual[3]]
    
    def check_52OL(numbers, letters):
        # both numbers incorrect, letters too early
        if numbers[0] in ['5', '2'] or numbers[1] in ['5', '2']:
            return False
        if letters[0] >= 'O' or letters[1] >= 'L':
            return False
        return True
    
    def check_20FP(numbers, letters):
        # one number correct in position, one incorrect, letters too early
        correct_pos = sum(1 for i in range(2) if numbers[i] == ['2', '0'][i])
        if correct_pos != 1:
            return False
        if letters[0] >= 'F' or letters[1] >= 'P':
            return False
        return True
    
    def check_08EX(numbers, letters):
        # one number correct wrong position, one too large
        if '8' in numbers:  # 8 is too large
            return False
        if numbers.count('0') != 1:  # 0 must be present but in wrong position
            return False
        return True
    
    def check_58RN(numbers, letters):
        # numbers too large, letters too early
        if int(numbers[0]) <= 5 or int(numbers[1]) <= 8:
            return False
        if letters[0] >= 'R' or letters[1] >= 'N':
            return False
        return True
    
    def check_64UM(numbers, letters):
        # one number wrong pos, one too large, one letter correct pos, one too early
        if '4' in numbers:  # 4 must be present but in wrong position
            pos_4 = numbers.index('4')
            if pos_4 == actual_numbers.index('4'):
                return False
        if 'U' not in letters or letters.index('U') != actual_letters.index('U'):
            return False
        return True
    
    def check_02PV(numbers, letters):
        # one number wrong pos, one incorrect
        if numbers.count('2') != 1:  # 2 must be present but in wrong position
            return False
        return True
    
    def check_62CO(numbers, letters):
        # both numbers incorrect, letters too early
        if '6' in numbers or '2' in numbers:
            return False
        if letters[0] >= 'C' or letters[1] >= 'O':
            return False
        return True
    
    def check_91CD(numbers, letters):
        # both numbers incorrect, letters too early
        if '9' in numbers or '1' in numbers:
            return False
        if letters[0] >= 'C' or letters[1] >= 'D':
            return False
        return True
    
    def check_96SR(numbers, letters):
        # numbers too large, one letter wrong pos, one too early
        if int(numbers[0]) <= 9 or int(numbers[1]) <= 6:
            return False
        if letters.count('S') != 1:
            return False
        return True

    checks = [
        check_52OL, check_20FP, check_08EX, check_58RN,
        check_64UM, check_02PV, check_62CO, check_91CD, check_96SR
    ]
    
    return all(check(numbers, letters) for check in checks)

# Generate all possible combinations
digits = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

possible_solutions = []

for d1, d2 in permutations(digits, 2):
    for l1, l2 in permutations(letters, 2):
        candidate = [d1, d2, l1, l2]
        if check_guess(candidate, candidate):
            possible_solutions.append(candidate)

print(possible_solutions)