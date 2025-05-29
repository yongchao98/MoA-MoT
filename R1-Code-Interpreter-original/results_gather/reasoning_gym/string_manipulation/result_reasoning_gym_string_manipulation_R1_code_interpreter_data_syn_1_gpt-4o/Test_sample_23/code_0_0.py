def transform_string(s):
    previous_states = set()
    
    while s not in previous_states:
        previous_states.add(s)
        original_s = s
        
        # Rule 1
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        # Rule 2
        elif s.startswith('cc'):
            s = s[2:]
        # Rule 3
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        # Rule 4
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Rule 5
        else:
            index = s.find('ca', 1)
            if index != -1:
                s = s[:index] + s[index+2:]
        
        # If no rule was applied, break the loop
        if s == original_s:
            break
    
    return s

# Initial string
initial_string = "aacbaca"
final_string = transform_string(initial_string)
print(final_string)