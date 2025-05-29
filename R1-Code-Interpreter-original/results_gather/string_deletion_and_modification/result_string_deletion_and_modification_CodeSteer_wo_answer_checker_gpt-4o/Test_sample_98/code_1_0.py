def transform_string(s):
    while True:
        original = s
        # Operation 1
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        # Operation 2
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # Operation 3
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        # Operation 4
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        # Operation 5
        elif 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]
        # Operation 6
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("bacccaccbca")
print(f'<<< "{final_string}" >>>')