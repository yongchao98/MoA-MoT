def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('cc'):
            s = s[2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        else:
            break
    return s

initial_string = "cccaababaaacaaaccb"
final_string = transform_string(initial_string)
print(final_string)