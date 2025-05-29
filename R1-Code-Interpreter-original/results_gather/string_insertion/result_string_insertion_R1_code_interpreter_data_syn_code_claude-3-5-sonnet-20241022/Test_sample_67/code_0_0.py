def apply_T_operation(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert to list for easier manipulation
    s = list(s)
    i = 0
    insertions = []  # Store (position, char) pairs
    
    # Scan the string
    while i <= len(s) - 4:
        substring = ''.join(s[i:i+4])
        if substring in patterns:
            # Store the insertion position and character
            insertions.append((i + 4, patterns[substring]))
        i += 1
    
    # Apply insertions from left to right
    # Need to account for previous insertions
    offset = 0
    for pos, char in insertions:
        s.insert(pos + offset, char)
        offset += 1
    
    result = ''.join(s)
    print(result)

# Test with the given string
s = "DCABCDCEABCAABCDCE"
apply_T_operation(s)