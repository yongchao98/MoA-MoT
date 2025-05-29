def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('cc'):
            s = s[2:]
        elif s.startswith('aa'):
            s = s[1:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        else:
            break
    return s

result = transform_string("ccbab")
print(result)