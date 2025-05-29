def apply_operations(s):
    while True:
        original = s
        # Operation 1
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        # Operation 2
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        # Operation 3
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        # Operation 4
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 5
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        # Operation 6
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

# Initial string
initial_string = "bcbababbcacbbabaaac"
final_string = apply_operations(initial_string)
print(final_string)