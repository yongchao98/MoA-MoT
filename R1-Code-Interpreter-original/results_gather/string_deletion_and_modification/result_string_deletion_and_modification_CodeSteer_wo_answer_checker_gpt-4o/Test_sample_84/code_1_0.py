def transform_string(s):
    while True:
        original = s
        # Operation 1
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 2
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        # Operation 3
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Operation 4
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        # Operation 5
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # Operation 6
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("abbbaaabcacbcbbcaba")
print(f'<<< "{final_string}" >>>')