def check_numbers(n1, n2):
    # Check all conditions for numbers
    nums = [n1, n2]
    
    # Condition 1: 98 both too large
    if n1 >= 9 or n2 >= 8: return False
    
    # Condition 2: 34 - one correct wrong pos, one too small
    c2_correct = (3 in nums) or (4 in nums)
    if not c2_correct: return False
    
    # Condition 3: 91 both incorrect
    if 9 in nums or 1 in nums: return False
    
    # Condition 4: 47 - one correct in pos, one too large
    if not ((n1 == 4 and n2 < 7) or (n2 == 7 and n1 < 7)): return False
    
    # Condition 5: 65 - one correct in pos, one too large
    if not ((n1 == 6 and n2 < 5) or (n2 == 5 and n1 < 6)): return False
    
    return True

def check_letters(l1, l2):
    # Check all conditions for letters
    letters = [l1, l2]
    
    # Condition 1: SO - one correct in pos, one incorrect
    so_count = sum(1 for x, y in zip(letters, ['S', 'O']) if x == y)
    if so_count != 1: return False
    
    # Condition 2: ZW - both too late
    if any(l >= 'Z' for l in letters): return False
    
    # Condition 3: HL - one correct in pos, one too early
    hl_count = sum(1 for x, y in zip(letters, ['H', 'L']) if x == y)
    if hl_count != 1: return False
    if any(l <= 'H' for l in letters if l != 'H'): return False
    
    # Condition 4: OB - both incorrect
    if 'O' in letters or 'B' in letters: return False
    
    # Condition 5: NA - both incorrect
    if 'N' in letters or 'A' in letters: return False
    
    return True

# Try all possible combinations
for n1 in range(10):
    for n2 in range(10):
        if n1 != n2 and check_numbers(n1, n2):
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 != l2 and check_letters(l1, l2):
                        print([str(n1), str(n2), l1, l2])