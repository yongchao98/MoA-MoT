def process_string(s):
    while True:
        initial = s
        # Operation 1
        if s.startswith('bb'):
            s = s[0] + s[2:]
        # Operation 2
        elif s.startswith('cc'):
            s = s[2:]
        # Operation 3
        elif s.endswith('bb'):
            s = s[:-2]
        # Operation 4
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        # Operation 5
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        # Operation 6
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        
        # If no operation was applied, break the loop
        if s == initial:
            break
    return s

final_string = process_string("cabbabbbabbbaabc")
print(final_string)