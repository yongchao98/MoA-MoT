def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('aa'):
            s = s[1:]
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        else:
            break
    return s

final_string = transform_string("bcbbabaa")
print(final_string)