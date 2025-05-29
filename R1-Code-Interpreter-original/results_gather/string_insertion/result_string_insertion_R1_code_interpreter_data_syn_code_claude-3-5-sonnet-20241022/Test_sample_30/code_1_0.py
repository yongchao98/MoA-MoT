def apply_T_operation(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    result = list(s)
    insertions = []
    
    # Scan and print found patterns
    for i in range(len(s)-3):
        substring = s[i:i+4]
        if substring in patterns:
            print(f"Found pattern {substring} at position {i}")
            insertions.append((i+4, patterns[substring]))
    
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    print(''.join(result))

apply_T_operation("ECDEAAEABCBDEAB")