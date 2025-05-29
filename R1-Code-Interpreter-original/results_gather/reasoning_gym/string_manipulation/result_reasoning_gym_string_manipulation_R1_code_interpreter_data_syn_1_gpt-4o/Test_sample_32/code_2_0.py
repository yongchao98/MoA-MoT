def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        else:
            break
    return s

initial_string = "accbccaabccaccaabbac"
final_string = transform_string(initial_string)
print(final_string)