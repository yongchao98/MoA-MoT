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
    i = 0
    insertions = []  # Store (position, char) pairs
    
    while i <= len(s) - 4:
        substring = s[i:i+4]
        if substring in patterns:
            # Store the insertion position and character
            insertions.append((i+4, patterns[substring]))
        i += 1
    
    # Apply insertions from left to right
    # Need to adjust positions for previous insertions
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    # Convert back to string
    final_result = ''.join(result)
    print(final_result)

# Test with the given string
s = "EBCDEEBCDEEDEABA"
apply_T_operation(s)