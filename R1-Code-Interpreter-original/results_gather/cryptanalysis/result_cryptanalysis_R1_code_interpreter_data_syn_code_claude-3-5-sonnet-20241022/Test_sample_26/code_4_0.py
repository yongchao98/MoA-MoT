def verify_clues(n1, n2, l1, l2):
    # Clue 1: 83VM
    # Both numbers wrong, one letter correct but wrong position
    if n1 in [8, 3] or n2 in [8, 3]:
        return False
    if not ('V' in [l1, l2] and 'V' != 'M'):  # V must be present
        return False

    # Clue 2: 40WL
    # Both numbers and letters wrong
    if n1 in [4, 0] or n2 in [4, 0] or 'W' in [l1, l2] or 'L' in [l1, l2]:
        return False

    # Clue 3: 50VD
    # Both numbers wrong, one letter correct but wrong position
    if n1 in [5, 0] or n2 in [5, 0] or 'D' in [l1, l2]:
        return False
    if 'V' not in [l1, l2]:  # V must be present
        return False

    # Clue 4: 40MC
    # Both numbers and letters wrong
    if n1 in [4, 0] or n2 in [4, 0] or 'M' in [l1, l2] or 'C' in [l1, l2]:
        return False

    # Clue 5: 18JC
    # One number correct and in position, one too large
    if n1 != 1:  # First number must be 1
        return False
    if 'J' in [l1, l2] or 'C' in [l1, l2]:
        return False

    # Clue 6: 41TW
    # One number correct but wrong position
    if n2 != 4:  # Second number must be 4 (since 1 is first)
        return False
    if 'T' in [l1, l2] or 'W' in [l1, l2]:
        return False

    # Clue 7: 07BF
    # Both numbers wrong, one letter correct but wrong position, one too early
    if n1 in [0, 7] or n2 in [0, 7] or 'B' in [l1, l2]:
        return False
    if l1 == 'F':  # F can't be first letter (too early)
        return False

    return True

# Test all valid letter combinations with fixed numbers [1,4]
valid_solutions = []
for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
        if l1 != l2:  # Letters can't repeat
            if verify_clues(1, 4, l1, l2):  # We know n1=1, n2=4 from clues 5 and 6
                valid_solutions.append([1, 4, l1, l2])

print(valid_solutions)