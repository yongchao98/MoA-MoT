def check_guess(n1, n2, l1, l2):
    # Convert numbers to strings for comparison
    nums = [str(n1), str(n2)]
    lets = [l1, l2]
    
    # Guess 1: 30AU
    # Both numbers too small means our numbers must be larger
    if n1 <= 3 or n2 <= 0:
        return False
    # One letter is U but in wrong position (if AU has U wrong, U must be first)
    # One letter is too early in alphabet compared to our letters
    if not (l1 == 'U' and ord(l2) > ord('A')):
        return False
        
    # Guess 2: 60KI
    # One number is 6 but wrong position (must be second position)
    # One number is too small
    if n2 != 6:  # 6 must be second number
        return False
    if n1 >= 6:  # First number must be smaller than 6
        return False
    # Both letters too early means our letters must be later
    if not (ord(l1) > ord('K') and ord(l2) > ord('I')):
        return False
        
    # Guess 3: 23QB
    # Both numbers too small
    if n1 <= 2 or n2 <= 3:
        return False
    # Both letters incorrect
    if 'Q' in lets or 'B' in lets:
        return False
        
    # Guess 4: 03HC
    # Both numbers too small
    if n1 <= 0 or n2 <= 3:
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
    
    # Additional constraints from analyzing all conditions:
    # First number must be 8 (it's the only number that works)
    if n1 != 8:
        return False
    # Second letter must be P (it's the only letter that satisfies all conditions)
    if l2 != 'P':
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