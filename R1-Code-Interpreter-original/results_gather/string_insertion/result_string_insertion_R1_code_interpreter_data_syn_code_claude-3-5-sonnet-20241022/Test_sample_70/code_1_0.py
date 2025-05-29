def apply_T_operation(s):
    # Convert to list for easier manipulation
    result = list(s)
    insertions = []  # Store (position, char) pairs
    
    # Scan the string
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        substring = s[i:i+4]
        if substring == 'ABCD':
            insertions.append((i+4, 'A'))
        elif substring == 'BCDE':
            insertions.append((i+4, 'B'))
        elif substring == 'CDEA':
            insertions.append((i+4, 'C'))
        elif substring == 'DEAB':
            insertions.append((i+4, 'D'))
        elif substring == 'EABC':
            insertions.append((i+4, 'E'))
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    print(''.join(result))

# Test the function
s = "DBABCDDABCDACDEAA"
apply_T_operation(s)