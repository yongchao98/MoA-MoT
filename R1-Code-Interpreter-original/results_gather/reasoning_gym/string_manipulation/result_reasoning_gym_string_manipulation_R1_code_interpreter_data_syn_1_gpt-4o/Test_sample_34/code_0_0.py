def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        else:
            break
    return s

initial_string = "bbabbbc"
final_string = transform_string(initial_string)
print(final_string)