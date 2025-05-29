def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if s.startswith('cc'):
            s = s[2:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        else:
            break
    return s

initial_string = "bbaaabccacb"
final_string = transform_string(initial_string)
print(final_string)