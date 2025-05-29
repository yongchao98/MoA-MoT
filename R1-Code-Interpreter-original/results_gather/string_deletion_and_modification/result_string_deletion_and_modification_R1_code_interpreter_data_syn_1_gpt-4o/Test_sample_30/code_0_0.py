def transform_string(s):
    while True:
        original = s
        # Operation 1
        if s.startswith('aa'):
            s = s[1:]
        # Operation 2
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 3
        if s.startswith('bb'):
            s = s[0] + s[2:]
        # Operation 4
        if s.startswith('cc'):
            s = s[2:]
        # Operation 5
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        # Operation 6
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

# Initial string
initial_string = "aaacbcaccacaacbcc"
# Transform the string
final_string = transform_string(initial_string)
print(final_string)