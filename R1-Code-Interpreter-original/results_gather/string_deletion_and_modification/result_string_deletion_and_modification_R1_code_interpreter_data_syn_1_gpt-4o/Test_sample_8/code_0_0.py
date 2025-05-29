def apply_operations(s):
    while True:
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        else:
            break
    return s

result = apply_operations("cbababccac")
print(result)