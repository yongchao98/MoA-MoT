def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.count('b') % 2 == 0 and 'b' in s:
            s += 'ab'
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif 'ca' in s[1:]:
            index = s.index('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.endswith('bb'):
            s = s[:-2]
        else:
            break
    return s

final_string = transform_string("acbaaaca")
print(final_string)