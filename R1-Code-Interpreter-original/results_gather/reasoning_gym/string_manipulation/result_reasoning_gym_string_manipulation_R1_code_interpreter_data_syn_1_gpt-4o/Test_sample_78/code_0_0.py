def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        else:
            break
    return s

initial_string = "cabbcbcccacbcbcaa"
final_string = transform_string(initial_string)
print(final_string)