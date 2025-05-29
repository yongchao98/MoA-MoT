def apply_T_operation(s):
    # Store the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Find all occurrences of all patterns
    insertions = []  # Will store tuples of (position, character)
    
    for i in range(len(s)):
        substr = s[i:i+4]  # Get 4 characters starting at position i
        if substr in patterns:
            insertions.append((i+4, patterns[substr]))
    
    # Sort insertions by position
    insertions.sort()
    
    # Apply insertions from left to right
    result = list(s)
    offset = 0  # Keep track of how many characters we've inserted
    
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    # Convert back to string and print
    print(''.join(result))

# Test with the given string
s = "EDCDEADBCDEEEABCBC"
apply_T_operation(s)