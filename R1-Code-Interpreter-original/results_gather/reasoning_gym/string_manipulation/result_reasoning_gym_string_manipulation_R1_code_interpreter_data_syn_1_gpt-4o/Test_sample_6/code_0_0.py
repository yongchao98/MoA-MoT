def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        else:
            break
    return s

initial_string = "bcaabcbbaccacbbbcb"
final_string = transform_string(initial_string)
print(final_string)