def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        else:
            break
    return s

final_string = transform_string("bcabbc")
print(final_string)