def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s

        # Rule 1
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            continue

        # Rule 2
        if s.startswith('bb'):
            s = s[0] + s[2:]
            continue

        # Rule 3
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            continue

        # Rule 4
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            continue

        # Rule 5
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            continue

        # Rule 6
        if s.count('b') % 2 == 0 and s.count('b') > 0:
            s = s + 'ab'
            continue

        # Rule 7
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
            continue

        # Rule 8
        if 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
            continue

        # If no rule was applied, break the loop
        if s == original_s:
            break

    return s

initial_string = "ccaabcccaabbcbaccaa"
final_string = transform_string(initial_string)
print(final_string)