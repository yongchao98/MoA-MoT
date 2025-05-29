def transform_string(s):
    while True:
        original = s
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
            print(f"Replaced 'acb' with 'bca': {s}")
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
            print(f"Replaced prefix 'ab' with 'ca': {s}")
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
            print(f"Replaced suffix 'ba' with 'ab': {s}")
        elif s.endswith('ca'):
            s = s[:-1]
            print(f"Removed last character from suffix 'ca': {s}")
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            print(f"Replaced prefix 'ca' with 'bb' and appended 'c': {s}")
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
            print(f"Removed 'ca' not at start: {s}")
        
        # If no changes were made, break the loop
        if s == original:
            break

    return s

final_string = transform_string("cacbcbcbac")
print(f'<<< "{final_string}" >>>')