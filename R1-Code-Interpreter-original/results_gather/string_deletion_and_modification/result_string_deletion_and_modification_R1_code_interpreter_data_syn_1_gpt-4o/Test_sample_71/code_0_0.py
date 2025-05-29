def apply_operations(s):
    while True:
        original = s
        # Operation 1
        if s.startswith('ab'):
            s = 'ca' + s[2:]
        # Operation 2
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        # Operation 3
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 4
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        # Operation 5
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Operation 6
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

final_string = apply_operations("bacccbbcacaaaacb")
print(final_string)