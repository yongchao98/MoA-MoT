def check_guess(guess, numbers, letters):
    # Check if a guess matches the given numbers and letters
    return (guess[0] == numbers[0] and 
            guess[1] == numbers[1] and 
            guess[2] == letters[0] and 
            guess[3] == letters[1])

def is_valid_solution(numbers, letters):
    # Convert numbers to strings for easier comparison
    guess = [str(numbers[0]), str(numbers[1]), letters[0], letters[1]]
    
    # Condition 1: 37ZI
    correct_pos = sum((a == b) for a, b in zip(['3','7','Z','I'], guess))
    correct_num = sum(x in guess[:2] for x in ['3','7'])
    if not (correct_pos == 1 and correct_num == 1): return False
    
    # Condition 2: 74VH
    correct_pos = sum((a == b) for a, b in zip(['7','4','V','H'], guess))
    correct_num = sum(x in guess[:2] for x in ['7','4'])
    if not (correct_pos == 0 and correct_num == 1): return False
    
    # Condition 3: 39SX
    if '3' in guess[:2] or '9' in guess[:2]: return False
    if not ((guess[2] == 'S' and guess[3] != 'X') or 
            (guess[3] == 'S' and guess[2] != 'X')): return False
    
    # Condition 4: 67GU
    correct_pos = sum((a == b) for a, b in zip(['6','7','G','U'], guess))
    correct_num = sum(x in guess[:2] for x in ['6','7'])
    if not (correct_pos == 1 and correct_num == 1): return False
    
    # Condition 5: 86AP
    if '8' in guess[:2] or '6' in guess[:2]: return False
    if not ((guess[2] == 'A' and guess[3] != 'P') or 
            (guess[3] == 'A' and guess[2] != 'P')): return False
    
    # Condition 6: 80LC
    if '8' in guess[:2] or '0' in guess[:2]: return False
    if 'L' in guess[2:] or 'C' in guess[2:]: return False
    
    # Condition 7: 38XZ
    if '3' in guess[:2] or '8' in guess[:2]: return False
    if 'Z' in guess[2:]: return False
    if not (('X' in guess[2:]) and (guess[2:].index('X') != ['X','Z'].index('X'))): return False
    
    # Condition 8: 96WO
    if '9' in guess[:2] or '6' in guess[:2]: return False
    if 'W' in guess[2:] or 'O' in guess[2:]: return False
    
    # Condition 9: 87BK
    if '8' in guess[:2]: return False
    correct_pos = sum((a == b) for a, b in zip(['8','7','B','K'], guess))
    if not (correct_pos == 1 and '7' in guess[:2]): return False
    
    # Condition 10: 84WL
    if '8' in guess[:2] or '4' in guess[:2]: return False
    if 'W' in guess[2:] or 'L' in guess[2:]: return False
    
    # Condition 11: 72NV
    if not (set(['7','2']) == set(guess[:2])): return False
    if 'N' in guess[2:] or 'V' in guess[2:]: return False
    
    return True

# Try all possible combinations
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2: continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2: continue
                if is_valid_solution([n1,n2], [l1,l2]):
                    print([str(n1), str(n2), l1, l2])