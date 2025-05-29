def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        else:
            break
    return s

initial_string = "cacabbccc"
final_string = transform_string(initial_string)
print(final_string)