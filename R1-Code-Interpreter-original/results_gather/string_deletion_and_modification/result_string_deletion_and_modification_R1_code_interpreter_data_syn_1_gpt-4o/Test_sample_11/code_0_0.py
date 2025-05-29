def apply_operations(s):
    while True:
        if s.startswith('aa'):
            s = s[1:]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            if index != -1:
                s = s[:index] + s[index+2:]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        else:
            break
    return s

final_string = apply_operations("abaabcaaca")
print(final_string)