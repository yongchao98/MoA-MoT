def check_guess(guess):
    n1, n2, l1, l2 = guess
    n1, n2 = int(n1), int(n2)
    
    # Check 18WV
    # For numbers: exactly one of (n1=1 or n2=8) must be true
    count_18 = (n1 == 1) + (n2 == 8)
    if count_18 != 1:
        return False
    
    # For letters in 18WV: exactly one of (l1=W or l2=V) must be true
    count_WV = (l1 == 'W') + (l2 == 'V')
    if count_WV != 1:
        return False
    
    # 76GH and 75LE: numbers must be larger than 7
    if n1 <= 7 and n1 != 1:  # exception for 1 if it's correct
        return False
    if n2 <= 7 and n2 != 8:  # exception for 8 if it's correct
        return False
    
    # 12BZ
    # For letters: exactly one of (l1=B or l2=Z) must be true
    count_BZ = (l1 == 'B') + (l2 == 'Z')
    if count_BZ != 1:
        return False
    
    # Additional letter constraints
    if l1 == l2:  # no repeating letters
        return False
    
    return True

solutions = []
# Try all possible combinations
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:  # no repeating numbers
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2:  # no repeating letters
                    continue
                guess = [str(n1), str(n2), l1, l2]
                if check_guess(guess):
                    solutions.append(guess)

print(solutions)