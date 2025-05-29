def apply_T_operation(s):
    # Define the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert to list for easier manipulation
    result = list(s)
    insertions = []  # Store (position, char) pairs
    
    # Scan the string
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        substring = s[i:i+4]
        if substring in patterns:
            insertions.append((i+4, patterns[substring]))
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    print(''.join(result))

# Test the string
test_string = "AEABCADEABAEABCBD"
apply_T_operation(test_string)