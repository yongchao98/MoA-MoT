def process_string(s):
    while True:
        original = s
        # Operation 1
        if s.endswith('bb'):
            s = s[:-2]
        # Operation 2
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        # Operation 3
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # Operation 4
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        # Operation 5
        elif s.endswith('ca'):
            s = s[:-1]
        # Operation 6
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

final_string = process_string("bbaaacccaaccabbca")
print(final_string)