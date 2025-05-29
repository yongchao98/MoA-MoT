def apply_operations(s):
    while True:
        original = s
        # Operation 1
        if s.startswith('bb'):
            s = s[0] + s[2:]
        
        # Operation 2
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        
        # Operation 3
        if s.endswith('bb'):
            s = s[:-2]
        
        # Operation 4
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        
        # Operation 5
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        
        # Operation 6
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

# Initial string
initial_string = "cbbaacbabccbcbba"
# Apply operations
final_string = apply_operations(initial_string)
print(final_string)