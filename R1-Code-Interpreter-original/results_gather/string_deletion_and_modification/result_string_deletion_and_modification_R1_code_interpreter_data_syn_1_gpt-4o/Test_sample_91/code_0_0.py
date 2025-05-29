def transform_string(s):
    while True:
        original = s
        if s.startswith('aa'):
            s = s[1:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('ca'):
            s = s[:-1]
        
        if s == original:
            break
    return s

result = transform_string("aabccbccbbccbb")
print(result)