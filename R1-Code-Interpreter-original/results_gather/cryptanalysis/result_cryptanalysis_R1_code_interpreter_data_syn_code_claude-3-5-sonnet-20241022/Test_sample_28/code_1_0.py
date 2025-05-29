def check_all_conditions(guess):
    # Basic conditions from previous code
    forbidden_nums = {1,5,8,0,3,4}
    forbidden_letters = {'Z','Q','I','K','U','M','J','Y','D'}
    
    n1, n2 = int(guess[0]), int(guess[1])
    l1, l2 = guess[2], guess[3]
    
    # Check basic conditions
    if n1 in forbidden_nums or n2 in forbidden_nums:
        return False
    if l1 in forbidden_letters or l2 in forbidden_letters:
        return False
    if n1 != 6:  # 6 must be in first position
        return False
    if 'R' not in [l1, l2]:  # R must be present
        return False
        
    # Additional conditions:
    # Test against guess 1 (15ZQ): both numbers and letters wrong
    if n1 == 1 or n2 == 5 or l1 == 'Z' or l2 == 'Q':
        return False
        
    # Test against guess 2 (51IK): letters too early in alphabet
    if ord(l1) <= ord('K') or ord(l2) <= ord('K'):
        return False
        
    # Test against guess 3 (21UM): one number correct but wrong position
    if n2 != 2:  # 2 must be in second position
        return False
        
    # Test against guess 4 (80JY)
    if n1 == 8 or n2 == 0 or l1 == 'J' or l2 == 'Y':
        return False
        
    # Test against guess 5 (59QY)
    if n1 == 5 or n2 == 9 or l1 == 'Q' or l2 == 'Y':
        return False
        
    # Test against guess 6 (65RI)
    # We already know 6 is correct in first position
    # R must be in wrong position from this guess
    if (l1 == 'R' and l2 == 'I') or (l1 == 'I' and l2 == 'R'):
        return False
        
    # Test against guess 7 (34MD)
    if n1 == 3 or n2 == 4 or l1 == 'M' or l2 == 'D':
        return False
        
    return True

# Generate all possible solutions
solutions = []
numbers = range(10)
letters = [chr(i) for i in range(65, 91)]  # A-Z

for n1 in numbers:
    for n2 in numbers:
        if n1 == n2:
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:
                    continue
                guess = f"{n1}{n2}{l1}{l2}"
                if check_all_conditions(guess):
                    solutions.append([str(n1), str(n2), l1, l2])

print(solutions)