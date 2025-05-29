def apply_T(s):
    result = list(s)
    i = 0
    insertions = 0  # Track number of insertions to adjust indices
    
    while i < len(s) - 3:
        substr = s[i:i+4]
        if substr == 'ABCD':
            result.insert(i + 4 + insertions, 'A')
            insertions += 1
        elif substr == 'BCDE':
            result.insert(i + 4 + insertions, 'B')
            insertions += 1
        elif substr == 'CDEA':
            result.insert(i + 4 + insertions, 'C')
            insertions += 1
        elif substr == 'DEAB':
            result.insert(i + 4 + insertions, 'D')
            insertions += 1
        elif substr == 'EABC':
            result.insert(i + 4 + insertions, 'E')
            insertions += 1
        i += 1
    
    print(''.join(result))

s = "EABCBDEABADEAB"
apply_T(s)