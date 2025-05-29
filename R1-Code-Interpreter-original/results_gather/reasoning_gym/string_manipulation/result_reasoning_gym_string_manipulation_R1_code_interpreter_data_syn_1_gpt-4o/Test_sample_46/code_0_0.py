def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('cc'):
            s = s[2:]
        else:
            break
    return s

initial_string = "abccbcbacac"
final_string = transform_string(initial_string)
print(final_string)