def transform_string(s):
    while True:
        original = s
        if s.startswith('aa'):
            s = s[1:]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('cc'):
            s = s[2:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        
        if s == original:
            break
    return s

result = transform_string("aaacbcaccacaacbcc")
print(result)