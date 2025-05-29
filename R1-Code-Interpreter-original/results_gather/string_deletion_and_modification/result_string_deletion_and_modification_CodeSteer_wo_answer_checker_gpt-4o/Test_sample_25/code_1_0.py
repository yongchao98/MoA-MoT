def transform_string(s):
    while True:
        original = s
        # Operation 1
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Operation 2
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        # Operation 3
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # Operation 4
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        # Operation 5
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 6
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("cbcccbcabcc")
print(f'<<< "{final_string}" >>>')