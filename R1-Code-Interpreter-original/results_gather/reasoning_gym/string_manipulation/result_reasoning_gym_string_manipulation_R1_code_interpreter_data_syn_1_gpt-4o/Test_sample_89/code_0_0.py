def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('ca'):
            s = s[:-1]
        else:
            break
    return s

initial_string = "acbcacc"
final_string = transform_string(initial_string)
print(final_string)