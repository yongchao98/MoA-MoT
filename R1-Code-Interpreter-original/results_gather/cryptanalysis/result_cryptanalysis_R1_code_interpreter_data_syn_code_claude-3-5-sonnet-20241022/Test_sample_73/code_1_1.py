def is_valid_password(pwd):
    n1, n2, l1, l2 = pwd
    
    # Must have 6 in first position and 8 in second position
    if n1 != '6' or n2 != '8':
        return False
    
    # Must have S in correct position (from condition 1)
    if l1 != 'S':
        return False
    
    # Must have Q in second letter position (from conditions)
    if l2 != 'Q':
        return False
    
    return True

numbers = ['6', '8']
letters = ['S', 'Q']

pwd = [numbers[0], numbers[1], letters[0], letters[1]]
print(pwd)