def apply_T_operation(s):
    # Define patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Keep track of insertions (position, character)
    insertions = []
    
    # Scan string from left to right
    for i in range(len(s)):
        # Check each pattern starting at current position
        for pattern, insert_char in patterns.items():
            if i + len(pattern) <= len(s):
                if s[i:i+len(pattern)] == pattern:
                    insertions.append((i + len(pattern), insert_char))
    
    # Sort insertions by position
    insertions.sort()
    
    # Apply insertions from left to right
    result = list(s)
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    return ''.join(result)

# Test with given string
s = "EAEABCBDEABACDEADB"
result = apply_T_operation(s)
print(result)