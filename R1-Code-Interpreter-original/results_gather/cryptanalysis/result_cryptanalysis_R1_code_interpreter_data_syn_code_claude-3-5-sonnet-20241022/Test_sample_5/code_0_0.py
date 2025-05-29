def check_numbers(n1, n2):
    # Check against all clues for numbers
    nums = [n1, n2]
    
    # Clue 1: 20 - both too small
    if n1 <= 2 or n2 <= 0:
        return False
    
    # Clue 2: 60 - one correct wrong pos, one too small
    c2_nums = {6, 0}
    correct_but_wrong = False
    too_small = False
    for n in nums:
        if n in c2_nums:
            correct_but_wrong = True
            c2_nums.remove(n)
    if not correct_but_wrong:
        return False
    
    # Clue 3: 39 - one correct in pos, one too large
    if n1 != 3 and n2 != 9:  # one must be 3 since it's correct in position
        return False
    
    # Clue 4: 09 - both incorrect
    if n1 == 0 or n1 == 9 or n2 == 0 or n2 == 9:
        return False
    
    return True

def check_letters(l1, l2):
    # Check against all clues for letters
    letters = {l1, l2}
    
    # Clue 1: KJ - both incorrect
    if 'K' in letters or 'J' in letters:
        return False
    
    # Clue 2: SC - one correct wrong pos, one incorrect
    if not ('S' in letters or 'C' in letters):
        return False
    
    # Clue 3: OQ - both incorrect
    if 'O' in letters or 'Q' in letters:
        return False
    
    # Clue 4: BV - one correct in pos, one too late
    if 'V' in letters:  # V is too late
        return False
    
    return True

# Try all possible combinations
numbers = range(10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
solution = None

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            if check_numbers(n1, n2):
                for l1 in letters:
                    for l2 in letters:
                        if l1 != l2:  # letters can't repeat
                            if check_letters(l1, l2):
                                solution = [str(n1), str(n2), l1, l2]

print(solution)