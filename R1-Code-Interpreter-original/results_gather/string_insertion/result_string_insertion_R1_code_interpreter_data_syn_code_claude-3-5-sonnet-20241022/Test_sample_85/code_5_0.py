def apply_T(s):
    # Convert string to list for easier manipulation
    result = list(s)
    length = len(s)
    i = 0
    insertions = 0
    
    # Process the original string without considering previous insertions
    matches = []
    while i < length - 3:
        substr = s[i:i+4]
        if substr == 'ABCD':
            matches.append((i, 'A'))
        elif substr == 'BCDE':
            matches.append((i, 'B'))
        elif substr == 'CDEA':
            matches.append((i, 'C'))
        elif substr == 'DEAB':
            matches.append((i, 'D'))
        elif substr == 'EABC':
            matches.append((i, 'E'))
        i += 1
    
    # Apply all insertions at once
    for pos, char in matches:
        result.insert(pos + 4 + insertions, char)
        insertions += 1
    
    print(''.join(result))

s = "EABCBDEABADEAB"
apply_T(s)