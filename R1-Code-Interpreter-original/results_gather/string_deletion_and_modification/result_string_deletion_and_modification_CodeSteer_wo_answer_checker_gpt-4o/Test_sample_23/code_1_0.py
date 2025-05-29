def transform_string(s):
    while True:
        original = s
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('cc'):
            s = s[2:]
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("babbcacacabba")
print(f'<<< "{final_string}" >>>')