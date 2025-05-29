def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        print(f"Current string: {s}")  # Debug: print the current state
        if s.startswith('bc'):
            s = s[2:] + 'aa'
            print(f"Applied Rule 1: {s}")  # Debug: print after applying Rule 1
        elif s.startswith('cc'):
            s = s[2:]
            print(f"Applied Rule 2: {s}")  # Debug: print after applying Rule 2
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
            print(f"Applied Rule 3: {s}")  # Debug: print after applying Rule 3
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            print(f"Applied Rule 4: {s}")  # Debug: print after applying Rule 4
        else:
            break
    return s

initial_string = "bbcaaaccb"
final_string = transform_string(initial_string)
print(f"Final string: {final_string}")