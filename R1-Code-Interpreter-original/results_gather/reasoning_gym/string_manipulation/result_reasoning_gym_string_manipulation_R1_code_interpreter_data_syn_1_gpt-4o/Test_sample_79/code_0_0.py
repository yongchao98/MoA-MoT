def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('aa'):
            s = s[1:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        else:
            break
    return s

initial_string = "bbabbcacacb"
final_string = transform_string(initial_string)
print(final_string)