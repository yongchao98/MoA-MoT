def transform_string(s):
    while True:
        original = s
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.count('b') % 2 == 0 and 'b' in s:
            s = s + 'ab'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("cbbcbbaabbcb")
print(f'<<< "{final_string}" >>>')