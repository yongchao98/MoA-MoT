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
        if s.endswith('ca'):
            s = s[:-1]
            continue

        # Rule 3
        if s.startswith('aa'):
            s = s[1:]
            continue

        # Rule 4
        if 'bca' in s:
            s = s.replace('bca', '', 1)
            continue

        # Rule 5
        if 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
            continue

        # If no rule was applied, break the loop
        if s == original_s:
            break

    return s

# Initial string
initial_string = "cccaaa"
final_string = transform_string(initial_string)
print(final_string)