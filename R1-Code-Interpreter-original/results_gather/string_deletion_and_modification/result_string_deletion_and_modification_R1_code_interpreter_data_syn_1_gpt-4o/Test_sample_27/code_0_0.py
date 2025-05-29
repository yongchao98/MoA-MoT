def transform_string(s):
    while True:
        initial_s = s
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.count('b') % 2 == 0 and 'b' in s:
            s = s + 'ab'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        
        # If no changes were made, break the loop
        if s == initial_s:
            break
    return s

# Initial string
initial_string = "cbbcbbaabbcb"
# Transform the string
final_string = transform_string(initial_string)
print(final_string)