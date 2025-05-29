def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('bb'):
            s = s[:-2]
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        else:
            break
    return s

initial_string = "cbbaccacabb"
final_string = transform_string(initial_string)
print(final_string)