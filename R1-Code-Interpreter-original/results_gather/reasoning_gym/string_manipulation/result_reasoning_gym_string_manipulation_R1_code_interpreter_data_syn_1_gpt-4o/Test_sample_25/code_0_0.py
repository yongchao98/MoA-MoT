def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        else:
            break
    return s

initial_string = "bccaac"
final_string = transform_string(initial_string)
print(final_string)