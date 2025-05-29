def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        else:
            break
    return s

final_string = transform_string("bcbaacbaa")
print(final_string)