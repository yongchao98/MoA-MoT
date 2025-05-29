def transform_string(s):
    while True:
        original = s
        # Operation 1
        if s.startswith('cc'):
            s = s[2:]
        # Operation 2
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Operation 3
        if s.startswith('aa'):
            s = s[1:]
        # Operation 4
        if s.endswith('ca'):
            s = s[:-1]
        # Operation 5
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # Operation 6
        if s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

result = transform_string("aaabcbbccbcaacc")
print(result)