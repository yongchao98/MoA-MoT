def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('aa'):
            s = s[1:]
        else:
            break
    return s

initial_string = "cbacbabccca"
final_string = transform_string(initial_string)
print(final_string)