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
    insertions = []  # Store (position, char) pairs for insertion
    
    # Scan the string
    for i in range(len(s)):
        # Check each pattern starting at position i
        for pattern, insert_char in patterns.items():
            if i + len(pattern) <= len(s):
                substring = s[i:i+len(pattern)]
                if substring == pattern:
                    insertions.append((i + len(pattern), insert_char))
    
    # Sort insertions by position (although they should already be in order)
    insertions.sort()
    
    # Apply insertions from right to left to maintain correct positions
    for pos, char in reversed(insertions):
        result.insert(pos, char)
    
    # Convert back to string
    result = ''.join(result)
    print(result)

# Test with the given string
s = "ABCDCCDEADCDEAEB"
apply_T_operation(s)