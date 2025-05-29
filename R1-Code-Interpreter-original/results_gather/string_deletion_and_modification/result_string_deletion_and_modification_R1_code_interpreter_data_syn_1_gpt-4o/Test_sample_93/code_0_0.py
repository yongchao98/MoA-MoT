def transform_string(s):
    while True:
        original = s
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        
        if s == original:
            break
    return s

final_string = transform_string("cbcbabcabcacbba")
print(final_string)