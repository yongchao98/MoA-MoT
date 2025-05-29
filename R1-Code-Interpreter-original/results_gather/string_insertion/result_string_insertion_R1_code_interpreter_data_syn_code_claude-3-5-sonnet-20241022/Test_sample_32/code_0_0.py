def apply_T_operation(s):
    # Define the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Find all occurrences of all patterns
    insertions = []
    for i in range(len(s)):
        substr = s[i:i+4]
        if substr in patterns:
            insertions.append((i+4, patterns[substr]))
    
    # Sort insertions by position
    insertions.sort()
    
    # Apply insertions from left to right
    result = list(s)
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    # Convert back to string
    result = ''.join(result)
    print(result)

# Test with the given string
s = "DBCDEBABCDDDEABAB"
apply_T_operation(s)