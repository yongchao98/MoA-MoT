def transform_string(s):
    while True:
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        else:
            break
    return s

initial_string = "abaacaabcab"
final_string = transform_string(initial_string)
print(final_string)