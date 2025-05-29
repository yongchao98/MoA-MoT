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
    
    # Scan the string once, from left to right
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        # Check the 4-character window at current position
        window = s[i:i+4]
        if window in patterns:
            insertions.append((i+4, patterns[window]))
    
    # Apply all insertions from left to right
    for offset, (pos, char) in enumerate(insertions):
        result.insert(pos + offset, char)
    
    print(''.join(result))

# Test with the given string
s = "CAEABCCEABCCBCDEEE"
apply_T_operation(s)