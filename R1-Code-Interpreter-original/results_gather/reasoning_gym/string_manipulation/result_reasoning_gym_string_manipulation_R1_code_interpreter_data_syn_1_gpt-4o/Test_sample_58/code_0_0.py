def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        else:
            break
    return s

initial_string = "cabcccbcbbc"
final_string = transform_string(initial_string)
print(final_string)