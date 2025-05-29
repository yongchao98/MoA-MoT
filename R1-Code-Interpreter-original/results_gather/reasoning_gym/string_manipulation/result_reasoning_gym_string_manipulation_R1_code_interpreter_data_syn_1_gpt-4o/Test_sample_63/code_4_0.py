def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        # Rule 1
        if len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        # Rule 2
        elif 'ca' in s[1:]:
            ca_index = s.index('ca', 1)
            s = s[:ca_index] + s[ca_index+2:]
        # Rule 3
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        else:
            break
    return s

initial_string = "baaacacba"
final_string = transform_string(initial_string)
print(final_string)