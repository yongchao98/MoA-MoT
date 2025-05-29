def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s

        # Rule 1
        if s.startswith('cc'):
            s = s[2:]
            continue

        # Rule 2
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            continue

        # Rule 3
        if s.endswith('bb'):
            s = s[:-2]
            continue

        # Rule 4
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            continue

        # Rule 5
        if 'bca' in s:
            s = s.replace('bca', '', 1)
            continue

        # If no rule was applied, break the loop
        if s == original_s:
            break

    return s

initial_string = "accccbabaabc"
final_string = transform_string(initial_string)
print(final_string)