def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s = s + 'ab'
        else:
            break
    return s

initial_string = "cacbcbbaa"
final_string = transform_string(initial_string)
print(final_string)