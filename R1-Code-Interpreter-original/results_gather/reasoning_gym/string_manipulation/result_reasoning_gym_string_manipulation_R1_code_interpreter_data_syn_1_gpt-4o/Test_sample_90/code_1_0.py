def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s  # Keep track of the original state before applying any rule
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('aa'):
            s = s[1:]
        
        # If no rule was applied, break the loop
        if s == original_s:
            break
    
    # Return the last unique state before any repetition
    return s

initial_string = "cbacbabccca"
final_string = transform_string(initial_string)
print(final_string)