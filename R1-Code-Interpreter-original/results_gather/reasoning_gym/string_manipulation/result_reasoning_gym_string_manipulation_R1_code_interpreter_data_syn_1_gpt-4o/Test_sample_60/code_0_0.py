def transform_string(s):
    seen_states = set()
    while True:
        if s in seen_states:
            break
        seen_states.add(s)
        
        if s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        else:
            break
    return s

initial_string = "cccccbcaccbaacccba"
final_string = transform_string(initial_string)
print(final_string)