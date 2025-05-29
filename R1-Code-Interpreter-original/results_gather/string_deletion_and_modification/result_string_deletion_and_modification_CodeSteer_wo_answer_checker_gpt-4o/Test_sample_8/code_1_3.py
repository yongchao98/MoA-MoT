def transform_string(s):
    while True:
        original = s
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
            print(f"Replaced 'aa' with 'cc': {s}")
        elif s.endswith('bb'):
            s = s[:-2]
            print(f"Deleted 'bb' from end: {s}")
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
            print(f"Replaced 'ab' with 'ca': {s}")
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
            print(f"Deleted 'bc' from start and appended 'aa': {s}")
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
            print(f"Replaced 'abc' with 'cab': {s}")
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
            print(f"Deleted 'bca': {s}")
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

final_string = transform_string("cbababccac")
print(f'<<< \"{final_string}\" >>>')