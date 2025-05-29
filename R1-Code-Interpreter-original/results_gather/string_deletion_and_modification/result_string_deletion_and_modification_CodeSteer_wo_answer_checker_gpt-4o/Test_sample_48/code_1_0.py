def transform_string(s):
    while True:
        original = s
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("bcbacacabbcaaacbcbc")
print(f'<<< "{final_string}" >>>')