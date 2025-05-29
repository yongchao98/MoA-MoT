def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        else:
            break
    return s

result = transform_string("bcbbbbabb")
print(result)