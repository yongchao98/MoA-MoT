def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('aa'):
            s = s[1:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        else:
            break
    return s

initial_string = "cccacaccacbbbcabaaa"
final_string = transform_string(initial_string)
print(final_string)