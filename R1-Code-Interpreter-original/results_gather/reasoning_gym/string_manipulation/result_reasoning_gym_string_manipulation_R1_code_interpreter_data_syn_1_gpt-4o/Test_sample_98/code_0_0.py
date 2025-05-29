def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.count('b') % 2 == 0 and 'b' in s:
            s = s + 'ab'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        else:
            break
    return s

initial_string = "bccbcbabacbabcccca"
final_string = transform_string(initial_string)
print(final_string)