def apply_operations(s):
    while True:
        original = s
        # Operation 1
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 2
        elif s.endswith('ca'):
            s = s[:-1]
        # Operation 3
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        # Operation 4
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Operation 5
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        # Operation 6
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        # If no operation was applied, break the loop
        if s == original:
            break
    return s

final_string = apply_operations("aaacaabcbabccba")
print(final_string)