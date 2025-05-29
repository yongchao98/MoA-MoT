def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        else:
            break
    return s

initial_string = "ccbaca"
final_string = transform_string(initial_string)
print(final_string)