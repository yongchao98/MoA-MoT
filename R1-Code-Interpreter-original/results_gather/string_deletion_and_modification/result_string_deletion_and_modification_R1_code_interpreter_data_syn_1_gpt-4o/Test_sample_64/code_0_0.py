def transform_string(s):
    while True:
        original = s
        # Operation 1
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        
        # Operation 2
        if len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        
        # Operation 3
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        
        # Operation 4
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
        
        # Operation 5
        if s.startswith('bb'):
            s = s[0] + s[2:]
        
        # Operation 6
        if s.startswith('aa'):
            s = s[1:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

initial_string = "baccabcaabcaabbcbca"
final_string = transform_string(initial_string)
print(final_string)