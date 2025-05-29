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
    for i in range(len(s)):
        # Check each pattern starting at position i
        for pattern, insert_char in patterns.items():
            if i + 4 <= len(s):  # Make sure we have enough characters left
                substring = s[i:i+4]
                if substring == pattern:
                    insertions.append((i+4, insert_char))
    
    # Sort insertions by position (although they should already be in order)
    insertions.sort()
    
    # Apply insertions from right to left to maintain correct positions
    for pos, char in reversed(insertions):
        result.insert(pos, char)
    
    return ''.join(result)

# Test with the given string
s = "EADEABDABCDCBCDE"
result = apply_T_operation(s)
print(result)