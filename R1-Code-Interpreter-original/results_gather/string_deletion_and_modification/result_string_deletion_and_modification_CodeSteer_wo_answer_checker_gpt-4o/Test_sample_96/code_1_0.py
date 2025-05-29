def transform_string(s):
    while True:
        original = s
        
        # Operation 1
        if s.startswith('bb'):
            s = s[0] + s[2:]
        
        # Operation 2
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        
        # Operation 3
        if 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        
        # Operation 4
        if s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        
        # Operation 5
        if s.startswith('ab'):
            s = 'ca' + s[2:]
        
        # Operation 6
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        
        # Check if any changes were made
        if s == original:
            break
    
    return s

final_string = transform_string("abacbbbcbbbacbcacba")
print(f'<<< "{final_string}" >>>')