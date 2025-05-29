def apply_operations(s):
    while True:
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        else:
            break
    return s

final_string = apply_operations("ccacacaacc")
print(final_string)