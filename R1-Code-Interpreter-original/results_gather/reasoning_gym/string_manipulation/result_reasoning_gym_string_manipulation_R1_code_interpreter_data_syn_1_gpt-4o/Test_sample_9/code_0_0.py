def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        else:
            break
    return s

initial_string = "baabccbcbbccaacab"
final_string = transform_string(initial_string)
print(final_string)