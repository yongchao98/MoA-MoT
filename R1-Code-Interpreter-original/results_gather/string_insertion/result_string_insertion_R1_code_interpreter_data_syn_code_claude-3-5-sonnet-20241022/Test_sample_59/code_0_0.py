def apply_T_operation(s):
    # Define patterns and their corresponding insertions
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
    
    # Check each position in the string
    for i in range(len(s)-3):
        window = s[i:i+4]
        if window in patterns:
            insertions.append((i+4, patterns[window]))
    
    # Apply insertions from left to right
    # Need to adjust positions after each insertion
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    # Convert back to string
    final_result = ''.join(result)
    print(final_result)

# Test with the given string
s = "ADEABEDEABBEABCBE"
apply_T_operation(s)