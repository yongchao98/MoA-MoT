from itertools import permutations

def check_guess(guess, actual):
    # Check if a guess matches all the feedback conditions
    return all([
        check_80YC(guess),
        check_58LT(guess),
        check_13QW(guess),
        check_36IL(guess),
        check_57CF(guess),
        check_46KW(guess),
        check_50GD(guess),
        check_05LX(guess),
        check_96SG(guess),
        check_45XI(guess),
        check_64DF(guess),
        check_31SE(guess),
        check_03LZ(guess)
    ])

def check_80YC(guess):
    return guess[0] != '8' and guess[1] != '0' and guess[2] != 'Y' and guess[3] != 'C'

def check_58LT(guess):
    nums = [guess[0], guess[1]]
    return ('5' in nums or '8' in nums) and '8' not in guess and 'L' not in guess and 'T' not in guess

def check_13QW(guess):
    return '1' not in guess and '3' not in guess and 'Q' not in guess and 'W' not in guess

def check_36IL(guess):
    return '3' not in guess and '6' not in guess and 'I' not in guess and 'L' not in guess

def check_57CF(guess):
    nums = [guess[0], guess[1]]
    return ('5' in nums or '7' in nums) and '7' not in guess and ('C' == guess[2] or 'C' == guess[3]) and 'F' not in guess

def check_46KW(guess):
    return '4' not in guess and '6' not in guess and 'K' not in guess and 'W' not in guess

def check_50GD(guess):
    nums = [guess[0], guess[1]]
    return ('5' in nums or '0' in nums) and guess[0] != '0' and 'G' not in guess and 'D' not in guess

def check_05LX(guess):
    return ('0' == guess[0] or '5' == guess[1]) and '0' not in [guess[1]] and 'L' not in guess and 'X' not in guess

def check_96SG(guess):
    return int(guess[0]) < 9 and int(guess[1]) < 9 and 'S' not in guess and 'G' not in guess

def check_45XI(guess):
    return ('4' == guess[0] or '5' == guess[1]) and '5' not in [guess[1]] and 'X' not in guess and 'I' not in guess

def check_64DF(guess):
    return '6' not in guess and '4' not in guess and ('D' == guess[2] or 'D' == guess[3]) and 'F' not in guess

def check_31SE(guess):
    return '3' not in guess and '1' not in guess and 'S' not in guess and 'E' not in guess

def check_03LZ(guess):
    return '0' not in guess and '3' not in guess and ('Z' == guess[2] or 'Z' == guess[3]) and 'L' not in guess

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Try all possible combinations
for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        guess = n1 + n2 + l1 + l2
                        if check_guess(guess, None):
                            print(f"Found solution: {[n1, n2, l1, l2]}")