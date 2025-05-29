def apply_operations(s):
    while True:
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            if index != -1:
                s = s[:index] + s[index+2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        else:
            break
    return s

final_string = apply_operations("bacccaccbca")
print(final_string)