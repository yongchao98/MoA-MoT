def check_numbers(n1, n2):
    # Convert numbers to strings for easier comparison
    num1, num2 = str(n1), str(n2)
    
    # Condition 1: 30AU - both numbers incorrect and too small
    if not (n1 > 3 and n2 > 0):
        return False
        
    # Condition 2: 60KI - one number correct but wrong position, one too small
    if not ((num1 == '0' and n2 < 6) or (num2 == '6' and n1 < 6) or 
            (num1 == '6' and n2 < 0) or (num2 == '0' and n1 < 6)):
        return False
        
    # Condition 3: 23QB - both numbers incorrect and too small
    if not (n1 > 2 and n2 > 3):
        return False
        
    # Condition 4: 03HC - both numbers incorrect and too small
    if not (n1 > 0 and n2 > 3):
        return False
        
    # Condition 5: 79MI - both numbers incorrect
    if num1 in '79' or num2 in '79':
        return False
        
    return True

def check_letters(l1, l2):
    # Condition 1: 30AU - one letter correct but wrong position, one too early
    if not ((l1 == 'U' and ord(l2) < ord('A')) or (l2 == 'A' and ord(l1) < ord('U'))):
        return False
        
    # Condition 2: 60KI - both letters too early
    if not (ord(l1) > ord('K') and ord(l2) > ord('I')):
        return False
        
    # Condition 3: 23QB - both letters incorrect
    if l1 in 'QB' or l2 in 'QB':
        return False
        
    # Condition 4: 03HC - both letters too early
    if not (ord(l1) > ord('H') and ord(l2) > ord('C')):
        return False
        
    # Condition 5: 79MI - both letters incorrect
    if l1 in 'MI' or l2 in 'MI':
        return False
        
    return True

def find_password():
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:  # Numbers can't repeat
                continue
            
            if not check_numbers(n1, n2):
                continue
                
            for l1 in [chr(i) for i in range(ord('A'), ord('Z')+1)]:
                for l2 in [chr(i) for i in range(ord('A'), ord('Z')+1)]:
                    if l1 == l2:  # Letters can't repeat
                        continue
                    
                    if check_letters(l1, l2):
                        return [str(n1), str(n2), l1, l2]
    return None

result = find_password()
print(result)