def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s

        # Rule 1: If the string starts with 'cc', remove the first two characters.
        if s.startswith('cc'):
            s = s[2:]
            continue

        # Rule 2: If the string suffix is 'cc', replace it with 'b' and prepend 'a' to the start.
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            continue

        # Rule 3: If the string prefix is 'aa', remove the first character.
        if s.startswith('aa'):
            s = s[1:]
            continue

        # Rule 4: If the string prefix is 'ca', replace it with 'bb' and append 'c' to the end.
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            continue

        # Rule 5: If the string starts with 'bb', remove the second character.
        if s.startswith('bb'):
            s = s[0] + s[2:]
            continue

        # If no rule was applied, break the loop
        if s == original_s:
            break

    return s

# Initial string
initial_string = "cccbcbbcabacbabcacbb"
final_string = transform_string(initial_string)
print(final_string)