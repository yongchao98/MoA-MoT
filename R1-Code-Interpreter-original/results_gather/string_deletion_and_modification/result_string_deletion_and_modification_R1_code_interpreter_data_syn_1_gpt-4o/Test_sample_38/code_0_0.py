def apply_operations(s):
    while True:
        original = s
        # Operation 1
        if s.endswith('ca'):
            s = s[:-1]
        # Operation 2
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        # Operation 3
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        # Operation 4
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        # Operation 5
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        # Operation 6
        elif s.startswith('aa'):
            s = s[1:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

final_string = apply_operations("baccaccaabcbabcacc")
print(final_string)