def transform_string(s):
    while True:
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        else:
            break
    return s

initial_string = "cacccaacaaaaa"
final_string = transform_string(initial_string)
print(final_string)