def verify_guess(guess, actual):
    n1, n2, l1, l2 = guess
    an1, an2, al1, al2 = actual
    
    # Count correct numbers in correct and wrong positions
    num_correct_pos = sum(g == a for g, a in zip([n1, n2], [an1, an2]))
    num_correct_wrong_pos = sum(g == a for g in [n1, n2] for a in [an1, an2]) - num_correct_pos
    
    # Count correct letters in correct and wrong positions
    let_correct_pos = sum(g == a for g, a in zip([l1, l2], [al1, al2]))
    let_correct_wrong_pos = sum(g == a for g in [l1, l2] for a in [al1, al2]) - let_correct_pos
    
    return (num_correct_pos, num_correct_wrong_pos, let_correct_pos, let_correct_wrong_pos)

def check_numbers_relation(n1, n2, test_n):
    return int(test_n) < min(int(n1), int(n2)) or int(test_n) > max(int(n1), int(n2))

def is_valid_solution(n1, n2, l1, l2):
    # Check condition 7: 56SM - both numbers correct but wrong positions
    if sorted([n1, n2]) != sorted(['5', '6']):
        return False
    
    # Check condition 10: 94AV - V must be in correct position
    if 'V' not in [l1, l2]:
        return False
        
    # Check condition 5: 79NP - N or P must be present
    if not ('N' in [l1, l2] or 'P' in [l1, l2]):
        return False
    
    # Verify specific guesses
    # Condition 1: 52FR
    g1 = verify_guess(['5', '2', 'F', 'R'], [n1, n2, l1, l2])
    if not (g1[0] + g1[1] == 1 and check_numbers_relation(n1, n2, '2')):
        return False
        
    # Condition 5: 79NP
    g5 = verify_guess(['7', '9', 'N', 'P'], [n1, n2, l1, l2])
    if not (g5[2] + g5[3] == 1):
        return False
        
    # Condition 8: 76BJ
    g8 = verify_guess(['7', '6', 'B', 'J'], [n1, n2, l1, l2])
    if not (g8[0] + g8[1] == 1):
        return False
        
    # Condition 10: 94AV
    g10 = verify_guess(['9', '4', 'A', 'V'], [n1, n2, l1, l2])
    if not (g10[2] == 1):
        return False
    
    return True

# Test all possible combinations
solutions = []
for n1, n2 in [('5','6'), ('6','5')]:
    for l1 in ['N', 'P', 'V']:
        for l2 in ['N', 'P', 'V']:
            if l1 != l2:
                if is_valid_solution(n1, n2, l1, l2):
                    solutions.append([n1, n2, l1, l2])

print("Solutions found:", solutions)