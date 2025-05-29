def check_guess(n1, n2, l1, l2):
    # Convert numbers to strings for comparison
    nums = [str(n1), str(n2)]
    lets = [l1, l2]
    
    # Guess 1: 30AU
    # Both numbers too small
    if not (n1 > 3 and n2 > 0):
        return False
    # One letter correct but wrong position, one too early
    if not ((l1 == 'U' or l2 == 'U') and ('A' not in [l1, l2])):
        return False
        
    # Guess 2: 60KI
    # One number is 6 but wrong position, other is too small
    if '6' not in nums or nums[0] == '6':  # 6 must be in second position
        return False
    # Both letters too early in alphabet
    if not (ord(l1) > ord('K') and ord(l2) > ord('I')):
        return False
        
    # Guess 3: 23QB
    # Both numbers too small
    if not (n1 > 2 and n2 > 3):
        return False
    # Both letters incorrect
    if 'Q' in lets or 'B' in lets:
        return False
        
    # Guess 4: 03HC
    # Both numbers too small
    if not (n1 > 0 and n2 > 3):
        return False
    # Both letters too early
    if not (ord(l1) > ord('H') and ord(l2) > ord('C')):
        return False
        
    # Guess 5: 79MI
    # Both numbers incorrect
    if '7' in nums or '9' in nums:
        return False
    # Both letters incorrect
    if 'M' in lets or 'I' in lets:
        return False
        
    return True

def find_password():
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:  # Numbers can't repeat
                continue
                
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 == l2:  # Letters can't repeat
                        continue
                    
                    if check_guess(n1, n2, l1, l2):
                        return [str(n1), str(n2), l1, l2]
    return None

result = find_password()
print(result)