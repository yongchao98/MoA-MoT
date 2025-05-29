def process_string(s):
    while True:
        original = s
        # Operation 1
        if 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]
        
        # Operation 2
        if s.endswith('bb'):
            s = s[:-2]
        
        # Operation 3
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        
        # Operation 4
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        
        # Operation 5
        if 'bca' in s:
            index = s.index('bca')
            s = s[:index] + s[index+3:]
        
        # Operation 6
        if s.startswith('cc'):
            s = s[2:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

# Initial string
initial_string = "acbbabbcabccacaccc"
# Process the string
final_string = process_string(initial_string)
print(final_string)