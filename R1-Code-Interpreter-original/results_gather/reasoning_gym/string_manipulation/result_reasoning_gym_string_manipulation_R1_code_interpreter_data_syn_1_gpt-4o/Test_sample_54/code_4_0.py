def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if len(s) > 1 and 'ca' in s[1:]:
            s = s[:s.index('ca', 1)] + s[s.index('ca', 1) + 2:]
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('aa'):
            s = s[1:]
        else:
            break
    return s

result = transform_string("cbcaa")
print(result)