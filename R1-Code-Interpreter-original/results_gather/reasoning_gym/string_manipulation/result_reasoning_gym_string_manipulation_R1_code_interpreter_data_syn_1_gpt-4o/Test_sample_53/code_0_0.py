def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        # Rule 1
        if s.count('b') % 2 == 0 and 'b' in s:
            s += 'ab'
            continue
        # Rule 2
        if s.startswith('ac'):
            s = 'zz' + s[2:]
            continue
        # Rule 3
        if s.startswith('aa'):
            s = s[1:]
            continue
        # Rule 4
        if s.startswith('bb'):
            s = s[0] + s[2:]
            continue
        # Rule 5
        if s.endswith('ca'):
            s = s[:-1]
            continue
        # Rule 6
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            continue
        # If no rule applies, break
        break
    return s

initial_string = "ccabaabbcac"
final_string = transform_string(initial_string)
print(final_string)