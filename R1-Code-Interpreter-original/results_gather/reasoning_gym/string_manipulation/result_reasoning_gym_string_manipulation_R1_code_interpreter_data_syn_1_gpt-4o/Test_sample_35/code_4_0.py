def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        else:
            break
    return s

final_string = transform_string("abbbbcb")
print(final_string)