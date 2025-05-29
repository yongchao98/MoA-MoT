def apply_operations(s):
    while True:
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        else:
            break
    return s

initial_string = "aaaaacabacbcbc"
final_string = apply_operations(initial_string)
print(final_string)