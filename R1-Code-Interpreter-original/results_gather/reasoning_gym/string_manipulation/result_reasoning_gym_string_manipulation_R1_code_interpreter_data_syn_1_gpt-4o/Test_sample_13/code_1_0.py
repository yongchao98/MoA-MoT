def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        else:
            break
    return s

initial_string = "bbbabbaabaaba"
final_string = transform_string(initial_string)
print(final_string)