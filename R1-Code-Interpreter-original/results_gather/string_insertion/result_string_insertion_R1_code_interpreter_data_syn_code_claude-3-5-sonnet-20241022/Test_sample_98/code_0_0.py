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
    
    # Scan string from left to right
    for i in range(len(s)):
        # Check each pattern
        for pattern, insert_char in patterns.items():
            # Check if pattern exists starting at position i
            if i + len(pattern) <= len(s):
                substring = s[i:i+len(pattern)]
                if substring == pattern:
                    insertions.append((i + len(pattern), insert_char))
    
    # Apply insertions from left to right
    # Sort by position to ensure correct order
    insertions.sort(key=lambda x: x[0])
    
    # Apply insertions (from left to right)
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    # Convert back to string
    result = ''.join(result)
    print(result)

# Test with given string
s = "EAABCDADEABECDEAAA"
apply_T_operation(s)