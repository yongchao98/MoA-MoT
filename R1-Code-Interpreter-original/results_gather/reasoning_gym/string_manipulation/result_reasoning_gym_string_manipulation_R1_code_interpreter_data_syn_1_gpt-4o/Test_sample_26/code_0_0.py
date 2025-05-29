def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s

        if s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('aa'):
            s = s[1:]

        if s == original_s:
            break

    return original_s

initial_string = "accaaabbbbb"
final_string = transform_string(initial_string)
print(final_string)