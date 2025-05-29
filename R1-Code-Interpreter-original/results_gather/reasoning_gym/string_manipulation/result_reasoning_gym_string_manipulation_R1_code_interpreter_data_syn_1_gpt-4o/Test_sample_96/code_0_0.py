def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        elif s.count('b') % 2 == 0 and 'b' in s:
            s += 'ab'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        else:
            break
    return s

initial_string = "bcbabaabccbc"
final_string = transform_string(initial_string)
print(final_string)