def check_all_conditions(guess):
    n1, n2 = int(guess[0]), int(guess[1])
    l1, l2 = guess[2], guess[3]
    
    # Basic number conditions
    if n1 != 6 or n2 != 2:  # We're certain about these numbers
        return False
    
    # Letter conditions
    forbidden_letters = {'Z','Q','I','K','U','M','J','Y','D'}
    if l1 in forbidden_letters or l2 in forbidden_letters:
        return False
    
    # Letters must be after K in alphabet
    if ord(l1) <= ord('K') or ord(l2) <= ord('K'):
        return False
    
    # From guess 6 (65RI):
    # R must be in a different position than in 65RI (where it was in position 3)
    # And one letter must be R
    if l1 != 'R' and l2 != 'R':
        return False
    if l1 == 'R' and l2 == 'I':  # Can't match the position in guess 6
        return False
    
    # Additional insight from guess 6:
    # Since we know 6 is correct and in position 1,
    # and R is correct but in wrong position (it was in position 3 in 65RI),
    # R must be in position 4 in our solution
    if l1 == 'R':  # R cannot be in position 3
        return False
    
    # The first letter (l1) must be after K but not any of the letters we've seen
    # in incorrect guesses
    seen_letters = {'Z','Q','I','K','U','M','J','Y','D','R'}
    if l1 in seen_letters:
        return False
    
    # From all guesses, we can deduce that the solution must have
    # letters that weren't tried in any incorrect positions
    tried_pairs = [
        ('Z','Q'), ('I','K'), ('U','M'), ('J','Y'), ('Q','Y'), ('R','I'), ('M','D')
    ]
    for pair in tried_pairs:
        if (l1, l2) == pair or (l2, l1) == pair:
            return False
    
    return True

# Generate all possible solutions
solutions = []
numbers = ['6', '2']  # We're certain about these
letters = [chr(i) for i in range(ord('L'), ord('Z')+1)]  # L-Z

for l1 in letters:
    for l2 in letters:
        if l1 == l2:
            continue
        guess = f"62{l1}{l2}"
        if check_all_conditions(guess):
            solutions.append([numbers[0], numbers[1], l1, l2])

print(solutions)