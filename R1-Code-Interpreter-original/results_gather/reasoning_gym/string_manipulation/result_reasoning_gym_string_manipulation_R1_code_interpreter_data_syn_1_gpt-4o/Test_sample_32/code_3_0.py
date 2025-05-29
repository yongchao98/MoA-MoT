def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        print(f"Current string: {s}")  # Debug statement
        seen_states.add(s)
        if s.startswith('bc'):
            s = s[2:] + 'aa'
            print("Applied Rule 1")  # Debug statement
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
            print("Applied Rule 2")  # Debug statement
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
            print("Applied Rule 3")  # Debug statement
        else:
            break
    return s

initial_string = "accbccaabccaccaabbac"
final_string = transform_string(initial_string)
print(final_string)