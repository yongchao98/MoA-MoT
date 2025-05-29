def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        else:
            break
    return s

initial_string = "acbccabacbcc"
final_string = transform_string(initial_string)
print(final_string)