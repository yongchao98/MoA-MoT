def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        else:
            break
    return s

initial_string = "ccccbbbabcab"
final_string = transform_string(initial_string)
print(final_string)