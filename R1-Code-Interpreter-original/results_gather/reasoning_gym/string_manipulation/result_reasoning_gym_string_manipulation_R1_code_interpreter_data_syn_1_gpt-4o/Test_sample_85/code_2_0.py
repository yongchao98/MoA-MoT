def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s  # Keep track of the original state before applying any rule
        if s.startswith('cc'):
            s = s[2:]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        
        # If no rule was applied, break the loop
        if s == original_s:
            break

    return s

initial_string = "aaabacbcaccc"
final_string = transform_string(initial_string)
print(final_string)