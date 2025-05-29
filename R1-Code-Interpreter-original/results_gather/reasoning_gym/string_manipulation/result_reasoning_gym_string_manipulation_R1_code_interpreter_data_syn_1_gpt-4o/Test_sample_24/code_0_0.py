def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('aa'):
            s = s[1:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        else:
            break
    return s

result = transform_string("bccbc")
print(result)