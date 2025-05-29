def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        else:
            break
    return s

initial_string = "ccbabaa"
final_string = transform_string(initial_string)
print(final_string)