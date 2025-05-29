def transform_string(s):
    while True:
        original = s
        # Operation 1
        if 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]
        
        # Operation 2
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        
        # Operation 3
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        
        # Operation 4
        if s.startswith('bb'):
            s = s[0] + s[2:]
        
        # Operation 5
        if 'bca' in s:
            index = s.index('bca')
            s = s[:index] + s[index+3:]
        
        # Operation 6
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("cabacaccbacabbcacca")
print(final_string)