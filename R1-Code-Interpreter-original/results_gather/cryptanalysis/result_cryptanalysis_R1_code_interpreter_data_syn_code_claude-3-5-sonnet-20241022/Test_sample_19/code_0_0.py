def check_numbers(n1, n2):
    # Check against all number conditions
    # 18WV: one correct position, one too small
    # 76GH: both too small
    # 75LE: both too small
    # 12BZ: both too small
    if n1 < 7 or n2 < 7:  # All numbers less than 7 are too small
        return False
    
    # Check 18WV condition
    if not ((n1 == 1 or n2 == 8) and (n1 != 1 or n2 != 8)):  # Only one of 1 or 8 is correct
        return False
    
    return True

def check_letters(l1, l2):
    # Check against all letter conditions
    # 18WV: one correct position, one too early
    # 76GH: both too early
    # 75LE: both too early
    # 12BZ: one correct position, one too early
    
    # Convert letters to numbers for comparison
    l1_num = ord(l1) - ord('A')
    l2_num = ord(l2) - ord('A')
    
    # Check if letters are after G, H, L, E (which are too early)
    if l1_num <= ord('H')-ord('A') or l2_num <= ord('H')-ord('A'):
        return False
    
    # Check 18WV condition
    wv_correct = (l1 == 'W' and 'V' > l2) or (l2 == 'V' and 'W' > l1)
    
    # Check 12BZ condition
    bz_correct = (l1 == 'Z' and 'B' > l2) or (l2 == 'Z' and 'B' > l1)
    
    return wv_correct or bz_correct

def find_password():
    # Try all possible combinations
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:  # Skip if numbers repeat
                continue
            if not check_numbers(n1, n2):
                continue
            
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 == l2:  # Skip if letters repeat
                        continue
                    if check_letters(l1, l2):
                        print([str(n1), str(n2), l1, l2])

find_password()