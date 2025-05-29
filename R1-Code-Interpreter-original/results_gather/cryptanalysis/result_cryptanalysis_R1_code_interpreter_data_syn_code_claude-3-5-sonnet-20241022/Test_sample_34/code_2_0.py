import string

def check_conditions(l1, l2):
    # Condition from 18VC
    # One letter correct and in position, one too early
    correct_pos = 0
    too_early = 0
    if l1 == 'V': correct_pos += 1
    if l2 == 'C': correct_pos += 1
    if l1 < 'V': too_early += 1
    if l2 < 'C': too_early += 1
    if not (correct_pos == 1 and too_early == 1):
        return False

    # Condition from 54OF
    # Both letters too early
    if not (l1 > 'O' and l2 > 'F'):
        return False

    # Condition from 36AV
    # V is correct but wrong position, A is too early
    has_v = (l1 == 'V' or l2 == 'V')
    if not has_v:
        return False
    if l1 == 'A' or l2 == 'A':
        return False

    # Condition from 42GE
    # Both letters too early
    if not (l1 > 'G' and l2 > 'E'):
        return False

    # Condition from 08NQ
    # Both letters incorrect
    if l1 in 'NQ' or l2 in 'NQ':
        return False

    # Additional check: letters must be different
    if l1 == l2:
        return False

    return True

valid_combinations = []
letters = string.ascii_uppercase

for l1 in letters:
    for l2 in letters:
        if check_conditions(l1, l2):
            valid_combinations.append(['3', '6', l1, l2])

print(valid_combinations)