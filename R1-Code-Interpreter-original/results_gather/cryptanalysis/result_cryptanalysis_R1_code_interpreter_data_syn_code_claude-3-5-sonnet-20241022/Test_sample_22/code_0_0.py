def check_guess(guess, numbers, letters):
    # Check if this combination satisfies all conditions
    
    # Condition 1: 83AJ
    correct_count = 0
    if 8 in numbers and numbers.index(8) != 0: correct_count += 1
    if 3 in numbers and numbers.index(3) != 1: correct_count += 1
    if correct_count != 1: return False
    if numbers[0] > 8 or numbers[1] > 8: return False
    if letters[0] <= 'A' or letters[1] <= 'J': return False
    
    # Condition 2: 72XZ
    if 7 in numbers or 2 in numbers: return False
    correct_count = 0
    if (letters[0] == 'X' and guess[2] == 'X') or (letters[1] == 'Z' and guess[3] == 'Z'):
        correct_count += 1
    if correct_count != 1: return False
    
    # Condition 3: 03LZ
    correct_count = 0
    if 0 in numbers and numbers.index(0) != 0: correct_count += 1
    if 3 in numbers and numbers.index(3) != 1: correct_count += 1
    if correct_count != 1: return False
    if numbers[0] < 0 or numbers[1] < 0: return False
    correct_count = 0
    if (letters[0] == 'L' and guess[2] == 'L') or (letters[1] == 'Z' and guess[3] == 'Z'):
        correct_count += 1
    if correct_count != 1: return False
    if letters[0] <= 'L' or letters[1] <= 'Z': return False
    
    # Condition 4: 25KR
    if 2 in numbers or 5 in numbers: return False
    if 'K' in letters or 'R' in letters: return False
    
    # Condition 5: 36TF
    correct_count = 0
    if (numbers[0] == 3 and guess[0] == '3') or (numbers[1] == 6 and guess[1] == '6'):
        correct_count += 1
    if correct_count != 1: return False
    if numbers[0] > 3 or numbers[1] > 6: return False
    if 'T' in letters or 'F' in letters: return False
    
    # Condition 6: 15JN
    if 1 in numbers or 5 in numbers: return False
    if letters[0] <= 'J' or letters[1] <= 'N': return False
    
    # Condition 7: 16FQ
    if 1 in numbers or 6 in numbers: return False
    correct_count = 0
    if ('F' in letters and letters.index('F') != guess[2:].index('F')): 
        correct_count += 1
    if correct_count != 1: return False
    if letters[0] <= 'F' or letters[1] <= 'Q': return False
    
    return True

# Try all possible combinations
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2: continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2: continue
                guess = f"{n1}{n2}{l1}{l2}"
                if check_guess(guess, [n1, n2], [l1, l2]):
                    print([str(n1), str(n2), l1, l2])