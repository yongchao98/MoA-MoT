def transform_string(s):
    while True:
        original = s
        if s.startswith('aa'):
            s = s[1:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('ca'):
            s = s[:-1]
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

# Initial string
initial_string = "aabccbccbbccbb"
# Transform the string
final_string = transform_string(initial_string)
print(final_string)