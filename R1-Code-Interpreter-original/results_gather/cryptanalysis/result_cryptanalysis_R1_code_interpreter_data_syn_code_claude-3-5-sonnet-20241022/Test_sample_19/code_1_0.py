def verify_guess(guess):
    # Check against all conditions
    
    # 18WV condition
    def check_18WV(n1, n2, l1, l2):
        nums_correct = (n1 == 1 or n2 == 8) and not (n1 == 1 and n2 == 8)  # one number correct
        nums_small = (n1 < 8 and n1 != 1) or (n2 < 8 and n2 != 8)  # one number too small
        letters_correct = (l1 == 'W' or l2 == 'V') and not (l1 == 'W' and l2 == 'V')  # one letter correct
        letters_early = (l1 < 'V' and l1 != 'W') or (l2 < 'W' and l2 != 'V')  # one letter too early
        return nums_correct and nums_small and letters_correct and letters_early
    
    # 76GH condition
    def check_76GH(n1, n2):
        return n1 > 7 and n2 > 7  # both numbers must be larger than 7
    
    # 12BZ condition
    def check_12BZ(l1, l2):
        letters_correct = (l1 == 'B' or l2 == 'Z') and not (l1 == 'B' and l2 == 'Z')
        letters_early = (l1 < 'Z' and l1 != 'B') or (l2 < 'B' and l2 != 'Z')
        return letters_correct and letters_early

    n1, n2, l1, l2 = guess
    n1, n2 = int(n1), int(n2)
    
    if not check_76GH(n1, n2):
        return False
        
    if not check_18WV(n1, n2, l1, l2):
        return False
        
    if not check_12BZ(l1, l2):
        return False
        
    return True

# Try all possible combinations
solutions = []
for n1 in range(8, 10):
    for n2 in range(8, 10):
        if n1 == n2:
            continue
        for l1 in 'VWXYZ':
            for l2 in 'VWXYZ':
                if l1 == l2:
                    continue
                guess = [str(n1), str(n2), l1, l2]
                if verify_guess(guess):
                    solutions.append(guess)

print(solutions)