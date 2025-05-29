def transform_string(s):
    while True:
        original = s
        # Operation 1: Remove 'ca' not at the start
        if 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]
        # Operation 2: Replace 'ac' at the start with 'zz'
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        # Operation 3: Replace 'ba' at the end with 'ab'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        # Operation 4: Replace 'cb' at the start with 'aa' and delete the last character
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        # Operation 5: Remove 'cc' at the start
        elif s.startswith('cc'):
            s = s[2:]
        # Operation 6: Remove the second character if the string starts with 'bb'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        
        # If no operation was applied, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("cbbcbcacccaabbcbcbca")
print(f'<<< "{final_string}" >>>')