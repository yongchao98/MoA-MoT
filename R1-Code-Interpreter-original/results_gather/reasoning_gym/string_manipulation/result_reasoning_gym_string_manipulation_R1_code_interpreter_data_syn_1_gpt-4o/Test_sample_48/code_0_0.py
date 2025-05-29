def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.count('b') % 2 == 0 and 'b' in s:
            s = s + 'ab'
        else:
            break
    return s

initial_string = "bccbbbcbbacba"
final_string = transform_string(initial_string)
print(final_string)