def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        else:
            break
    return s

initial_string = "bbcaaaccb"
final_string = transform_string(initial_string)
print(final_string)