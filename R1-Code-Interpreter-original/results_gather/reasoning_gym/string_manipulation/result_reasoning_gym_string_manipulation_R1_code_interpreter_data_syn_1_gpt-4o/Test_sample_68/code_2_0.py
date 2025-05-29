def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        print(f"Current string: {s}")
        if s.count('b') % 2 == 0 and 'b' in s:
            s += 'ab'
            print(f"Applied Rule 1: {s}")
        elif s.endswith('bb'):
            s = s[:-2]
            print(f"Applied Rule 2: {s}")
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
            print(f"Applied Rule 3: {s}")
        else:
            break
    return s

initial_string = "bbacbcbbc"
final_string = transform_string(initial_string)
print(f"Final string: {final_string}")