def apply_T_operation(s):
    result = s
    i = 0
    insertions = []
    
    while i <= len(s) - 4:
        substr = s[i:i+4]
        if substr == 'ABCD':
            insertions.append((i+4, 'A'))
        elif substr == 'BCDE':
            insertions.append((i+4, 'B'))
        elif substr == 'CDEA':
            insertions.append((i+4, 'C'))
        elif substr == 'DEAB':
            insertions.append((i+4, 'D'))
        elif substr == 'EABC':
            insertions.append((i+4, 'E'))
        i += 1
    
    # Apply insertions from left to right
    result = list(s)
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
        
    print(''.join(result))

s = "BCDEEBCDECABCDDD"
apply_T_operation(s)