def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('cc'):
            s = s[2:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        else:
            break
    return s

initial_string = "cbbaccccbcccbacbabb"
final_string = transform_string(initial_string)
print(final_string)