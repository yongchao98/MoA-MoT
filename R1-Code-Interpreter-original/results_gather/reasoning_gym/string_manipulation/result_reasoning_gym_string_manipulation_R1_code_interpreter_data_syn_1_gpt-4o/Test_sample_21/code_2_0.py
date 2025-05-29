def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        else:
            break
    return s

initial_string = "bccbccbcbca"
final_string = transform_string(initial_string)
print(final_string)