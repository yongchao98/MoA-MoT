def verify_guess(guess, actual):
    n1, n2, l1, l2 = guess
    an1, an2, al1, al2 = actual
    
    # Numbers analysis
    nums_correct_pos = sum(g == a for g, a in zip([n1, n2], [an1, an2]))
    nums_wrong_pos = sum(g == a for g in [n1, n2] for a in [an1, an2]) - nums_correct_pos
    nums_too_small = sum(1 for n in [n1, n2] if int(n) < min(int(an1), int(an2)))
    nums_too_large = sum(1 for n in [n1, n2] if int(n) > max(int(an1), int(an2)))
    
    # Letters analysis
    lets_correct_pos = sum(g == a for g, a in zip([l1, l2], [al1, al2]))
    lets_wrong_pos = sum(g == a for g in [l1, l2] for a in [al1, al2]) - lets_correct_pos
    
    return (nums_correct_pos, nums_wrong_pos, nums_too_small, nums_too_large,
            lets_correct_pos, lets_wrong_pos)

def is_valid_solution(n1, n2, l1, l2):
    # Condition 1: 52FR
    # One number correct but wrong position, one too small, both letters wrong
    c1 = verify_guess(['5', '2', 'F', 'R'], [n1, n2, l1, l2])
    if not (c1[1] == 1 and c1[2] == 1 and c1[4] == 0 and c1[5] == 0):
        return False

    # Condition 2: 02IT
    # Both numbers too small, both letters wrong
    c2 = verify_guess(['0', '2', 'I', 'T'], [n1, n2, l1, l2])
    if not (c2[2] == 2 and c2[4] == 0 and c2[5] == 0):
        return False

    # Condition 5: 79NP
    # Both numbers too large, one letter correct but wrong position
    c5 = verify_guess(['7', '9', 'N', 'P'], [n1, n2, l1, l2])
    if not (c5[3] == 2 and c5[5] == 1):
        return False

    # Condition 7: 56SM
    # Both numbers correct but wrong positions
    c7 = verify_guess(['5', '6', 'S', 'M'], [n1, n2, l1, l2])
    if not (c7[1] == 2 and c7[4] == 0 and c7[5] == 0):
        return False

    # Condition 8: 76BJ
    # One number correct but wrong position, one too large
    c8 = verify_guess(['7', '6', 'B', 'J'], [n1, n2, l1, l2])
    if not (c8[1] == 1 and c8[3] == 1):
        return False

    # Condition 10: 94AV
    # Both numbers wrong, one letter correct and in correct position
    c10 = verify_guess(['9', '4', 'A', 'V'], [n1, n2, l1, l2])
    if not (c10[0] == 0 and c10[1] == 0 and c10[4] == 1 and c10[5] == 0):
        return False

    # Additional conditions:
    # V must be in correct position (from condition 10)
    if l2 != 'V':  # V must be in second letter position
        return False

    # N must be in first letter position (from condition 5)
    if l1 != 'N':
        return False

    return True

# Test all possible combinations
solutions = []
for n1, n2 in [('5','6'), ('6','5')]:
    for l1 in ['N', 'P']:
        if is_valid_solution(n1, n2, l1, 'V'):
            solutions.append([n1, n2, l1, 'V'])

print("Solutions found:", solutions)