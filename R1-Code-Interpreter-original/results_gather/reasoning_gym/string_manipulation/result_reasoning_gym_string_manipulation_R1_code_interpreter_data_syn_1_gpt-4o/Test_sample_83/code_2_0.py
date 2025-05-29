def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        # Rule 1
        if s.count('b') % 2 == 0 and 'b' in s:
            s += 'ab'
            continue
        # Rule 2
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
            continue
        # Rule 3
        if len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            continue
        # Rule 4
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            continue
        # Rule 5
        if 'bca' in s:
            s = s.replace('bca', '', 1)
            continue
        # Rule 6
        if s.endswith('ca'):
            s = s[:-1]
            continue
        # Rule 7
        if s.startswith('bb'):
            s = s[0] + s[2:]
            continue
        break
    return s

initial_string = "aacabcbccbcbba"
final_string = transform_string(initial_string)
print(final_string)