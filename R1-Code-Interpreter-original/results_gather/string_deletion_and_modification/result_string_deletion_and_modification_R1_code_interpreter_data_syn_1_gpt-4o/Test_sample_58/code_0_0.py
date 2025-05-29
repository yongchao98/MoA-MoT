def transform_string(s):
    while True:
        original = s
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        
        if s == original:
            break
    return s

final_string = transform_string("cacbcbcbac")
print(final_string)