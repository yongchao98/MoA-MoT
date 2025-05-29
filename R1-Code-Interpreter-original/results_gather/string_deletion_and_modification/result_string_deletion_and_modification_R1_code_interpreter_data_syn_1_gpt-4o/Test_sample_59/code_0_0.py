def transform_string(s):
    while True:
        original = s
        # Operation 1
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        # Operation 2
        elif 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]
        # Operation 3
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        # Operation 4
        elif s.endswith('ca'):
            s = s[:-1]
        # Operation 5
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        # Operation 6
        elif s.startswith('aa'):
            s = s[1:]
        # If no changes were made, break the loop
        if s == original:
            break
    return s

result = transform_string("abcbabaccacaba")
print(result)