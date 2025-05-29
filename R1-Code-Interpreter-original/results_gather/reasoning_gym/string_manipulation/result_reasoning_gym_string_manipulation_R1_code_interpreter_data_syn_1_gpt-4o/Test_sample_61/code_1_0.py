def transform_string(s):
    seen_states = []
    
    while True:
        if s in seen_states:
            return seen_states[-1]  # Return the state before the repeated state
        seen_states.append(s)
        
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('aa'):
            s = s[1:]
        else:
            break
    
    return s

initial_string = "ccbaabaabaacaccaccbb"
final_string = transform_string(initial_string)
print(final_string)