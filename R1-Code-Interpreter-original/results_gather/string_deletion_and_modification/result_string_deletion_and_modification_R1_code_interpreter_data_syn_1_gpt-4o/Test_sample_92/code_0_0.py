def apply_operations(s):
    while True:
        original = s
        # Operation 1
        if 'ca' in s[1:]:
            s = s.replace('ca', '', 1)
        # Operation 2
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Operation 3
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        # Operation 4
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        # Operation 5
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        # Operation 6
        elif s.startswith('cc'):
            s = s[2:]
        # If no operation was applied, break the loop
        if s == original:
            break
    return s

final_string = apply_operations("acbbcabcccbcc")
print(final_string)