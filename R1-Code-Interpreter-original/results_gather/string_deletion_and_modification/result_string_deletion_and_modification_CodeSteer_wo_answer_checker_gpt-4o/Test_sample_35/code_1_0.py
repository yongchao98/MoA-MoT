def transform_string(s):
    while True:
        original = s

        # 1. Remove 'ca' (not at the start)
        if 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]

        # 2. Append 'ab' if even number of 'b's
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'

        # 3. Remove 'bca'
        elif 'bca' in s:
            index = s.index('bca')
            s = s[:index] + s[index+3:]

        # 4. Replace 'cb' prefix
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]

        # 5. Replace 'ab' prefix
        elif s.startswith('ab'):
            s = 'ca' + s[2:]

        # 6. Replace 'cc' suffix
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'

        # If no changes were made, break the loop
        if s == original:
            break

    return s

final_string = transform_string("acaccbcbaacccbcccaac")
print(f'<<< "{final_string}" >>>')