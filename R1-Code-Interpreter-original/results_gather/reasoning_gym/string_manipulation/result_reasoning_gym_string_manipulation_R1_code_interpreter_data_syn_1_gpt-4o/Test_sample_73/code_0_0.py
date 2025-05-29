def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if s.endswith('ca'):
            s = s[:-1]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        else:
            break
    return s

initial_string = "caabccbbabaaacbb"
final_string = transform_string(initial_string)
print(final_string)