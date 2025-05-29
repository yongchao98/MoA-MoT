def transform_string(s):
    while True:
        original = s
        # Operation 1
        if 'ca' in s[1:]:
            s = s[:s.index('ca', 1)] + s[s.index('ca', 1) + 2:]
        # Operation 2
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 3
        elif s.endswith('bb'):
            s = s[:-2]
        # Operation 4
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # Operation 5
        elif s.startswith('aa'):
            s = s[1:]
        # Operation 6
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        
        # If no operations were applied, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("aacbbacacabbba")
print(f'<<< "{final_string}" >>>')