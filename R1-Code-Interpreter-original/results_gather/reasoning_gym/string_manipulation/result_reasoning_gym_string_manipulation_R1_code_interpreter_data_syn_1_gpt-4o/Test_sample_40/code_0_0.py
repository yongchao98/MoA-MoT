def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('aa'):
            s = s[1:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        else:
            break
    return s

initial_string = "ccabccc"
final_string = transform_string(initial_string)
print(final_string)