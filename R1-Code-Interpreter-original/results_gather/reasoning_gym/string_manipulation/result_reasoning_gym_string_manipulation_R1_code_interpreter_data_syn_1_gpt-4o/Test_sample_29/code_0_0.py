def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        else:
            break
    return s

initial_string = "acabbbcabbbabbbabcac"
final_string = transform_string(initial_string)
print(final_string)