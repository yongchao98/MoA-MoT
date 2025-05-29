def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('cc'):
            s = s[2:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif 'ca' in s[1:]:
            ca_index = s.find('ca', 1)
            s = s[:ca_index] + s[ca_index+2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        else:
            break
    return s

initial_string = "cabcaccaaccc"
final_string = transform_string(initial_string)
print(final_string)