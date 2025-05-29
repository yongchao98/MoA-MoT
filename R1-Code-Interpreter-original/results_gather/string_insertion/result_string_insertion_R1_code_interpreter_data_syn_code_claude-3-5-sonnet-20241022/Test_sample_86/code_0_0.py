def apply_T_operation(s):
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
            if i + 4 <= len(s):  # Make sure we have enough characters to check
                substring = s[i:i+4]
                if substring == pattern:
                    insertions.append((i+4, insert_char))
    
    # Apply insertions from left to right
    # Sort by position to ensure correct order
    insertions.sort(key=lambda x: x[0])
    
    # Apply insertions (need to account for previous insertions)
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    final_string = ''.join(result)
    print(final_string)

# Test with the given string
s = "EABCBABCDADEABDB"
apply_T_operation(s)