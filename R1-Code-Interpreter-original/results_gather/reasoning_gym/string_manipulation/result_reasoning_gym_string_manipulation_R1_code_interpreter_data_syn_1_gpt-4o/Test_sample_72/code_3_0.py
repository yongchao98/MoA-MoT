def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        else:
            break
    return s

initial_string = "acabb"
final_string = transform_string(initial_string)
print(final_string)