def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.count('b') % 2 == 0 and 'b' in s:
            s += 'ab'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        else:
            break
    return s

initial_string = "bbacbcbbc"
final_string = transform_string(initial_string)
print(final_string)