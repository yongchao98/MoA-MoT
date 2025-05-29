def transform_string(s):
    while True:
        original = s
        if s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("bbcaaccbbccbbbaa")
print(f'<<< "{final_string}" >>>')