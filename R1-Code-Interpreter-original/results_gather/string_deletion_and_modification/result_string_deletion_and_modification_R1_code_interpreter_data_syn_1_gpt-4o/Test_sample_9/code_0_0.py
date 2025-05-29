def transform_string(s):
    while True:
        original = s
        # Operation 1
        if 'ca' in s[1:]:
            s = s[:s.index('ca', 1)] + s[s.index('ca', 1) + 2:]
        # Operation 2
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        # Operation 3
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        # Operation 4
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        # Operation 5
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # Operation 6
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        # If no operation was applied, break the loop
        if s == original:
            break
    return s

final_string = transform_string("bbcbccacaacbbab")
print(final_string)